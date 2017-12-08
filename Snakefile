#!/usr/bin/env python

configfile: "config.yaml"
SAMPLES = config["samples"]
IP = [v["ip"] for k,v in SAMPLES.items()]

def uniq(input):
    output=[]
    for x in input:
        if x not in output:
            output.append(x)
    return output



localrules: all,
   	    make_barcode_file


rule all:
	input:
		"barcodes.tsv",
		"qual_ctrl/raw",
		expand("qual_ctrl/trimmed-{sample}", sample=SAMPLES),
		expand("alignment/fastq/unaligned-{sample}.fastq.gz", sample=SAMPLES),
		expand("alignment/fastq/aligned-{sample}.fastq.gz", sample=SAMPLES),
		expand("fastq/trimmed/{sample}.trimmed.fastq.gz", sample=SAMPLES),
		expand("alignment/{sample}.bam.bai", sample = SAMPLES),
		expand("cross_correlation/fragment_length-{sample}.txt", sample = SAMPLES),
		"alignment/bowtie_read_number_summary.txt",
		expand("wigfiles/{sample}_normSI.wig", sample = SAMPLES),
		expand("wigfiles/Spom/{sample}_Spom_normSI.wig", sample = SAMPLES),
		"reads_stats/spikein_proportion_all.txt",
		expand("reads_stats/plots/spike_in_proportion_{ip}.{type}", ip=uniq(IP), type = config["imagetype"]["plot_spikein_percentage"])

		
		
#Make a tab separated file with column 1 containing name of lib and column 2 containing the barcode sequence
rule make_barcode_file:
	output:
		"barcodes.tsv"
	params:
		bc = config["samples"]
	run:
		with open(output[0],"w") as f:
			for x in params.bc:
				f.write('\t'.join((x,params.bc[x]["barcode"])) + '\n')
#fastqc analysis				
rule fastqc_raw:
	input:
		config["fastq"]
	output:
		"qual_ctrl/raw"
	threads: config["threads"]
	log: "logs/fastqc/fastqc-raw.log"
	shell: """
	mkdir {output}
	(fastqc -o {output} --noextract -t {threads} {input}) &> {log}
	"""
#Demultiplex fastq file, allow one mismatch in barcode as long as they are unique
rule demultiplex:
	input:
		fastq = config["fastq"],
		barcode = "barcodes.tsv"
	output:
		expand("fastq/{sample}.fastq.gz", sample=SAMPLES)
	log: "logs/demultiplex.log"
	params:
		mismatch = config["demultiplex"]["mismatch"]
	shell:"""
	(fastq-multx -B {input.barcode} -b {input.fastq} -m {params.mismatch} -o fastq/%.fastq.gz) &> {log}
	"""
#Cutadapt - trim 1 base from 5' end of every file
#Also trims low quality reads from the 3' end and removes reads that are less than 'n' nucleotides in length after trimming
rule cutadapt:
	input:
		"fastq/{sample}.fastq.gz"
	output:
		temp("fastq/trimmed/{sample}.trimmed.fastq")
	params:
		qual_cutoff = config["cutadapt"]["qual_cutoff"],
		min_length = config["cutadapt"]["min_length"]
	log: "logs/cutadapt/cutadapt-{sample}.log"
	shell:"""
	(cutadapt -u 1 --nextseq-trim={params.qual_cutoff} -m {params.min_length} -o {output} {input}) &> {log}
	"""
#fastQC on demultiplexed and trimmed reads
rule fastqc_processed:
	input:
		"fastq/trimmed/{sample}.trimmed.fastq"
	output:
		"qual_ctrl/trimmed-{sample}"
	threads: config["threads"]
	log: "logs/fastqc/fastqc-trimmed-{sample}.log"
	shell:"""
	mkdir {output}
	(fastqc -o {output} --noextract -t {threads} {input}) &> {log}
	"""
#Alignment to the genome using bowtie
basename = "ScSp"
rule bowtie:
	input:
		fastq = "fastq/trimmed/{sample}.trimmed.fastq"
	output:
		sam = temp("alignment/{sample}.sam"),
		unaligned = temp("alignment/fastq/unaligned-{sample}.fastq"),
		aligned = temp("alignment/fastq/aligned-{sample}.fastq")
	params:
		index = "genome/" + basename
	threads: config["threads"]
	log: "logs/bowtie/bowtie-align-{sample}.log"
	shell:"""
	(bowtie2 -p {threads} -x {params.index} -U {input.fastq} --al {output.aligned} --un {output.unaligned} -S {output.sam}) &> {log} 
	"""
#Convert SAM file to BAM and sort the file
rule samtools_sort:
	input:
		"alignment/{sample}.sam"
	output:
		"alignment/{sample}.bam"
	log: "logs/samtools_sort/samtools-sort-{sample}.log"
	threads: config["threads"]
	shell:"""
	(samtools view -buh -q 3 {input} | samtools sort -T {wildcards.sample} -@ {threads} -o {output} -) &> {log}
	"""
#Gather summary statistics for alignment
rule bowtie_summary:
	input:
		expand("alignment/{sample}.bam", sample=SAMPLES)
	output:
		number = "alignment/bowtie_read_number_summary.txt",
		percentage = "alignment/bowtie_read_percentage_summary.txt"
	script: "scripts/extract_bowtie_log.R"
#Samtools indexing - make bam.bai files
rule samtools_index:
	input:
		"alignment/{sample}.bam"
	output:
		"alignment/{sample}.bam.bai"
	log: "logs/samtools_index/samtools_index-{sample}.log"
	shell:"""
	(samtools index -b {input}) &> {log}
	"""
#Zip the trimmed fastq and unaligned/aligned fastq files
rule gzip_loose_fastq:
	input:
		"alignment/{sample}.sam",
		trim = "fastq/trimmed/{sample}.trimmed.fastq",
		aligned = "alignment/fastq/unaligned-{sample}.fastq",
		unaligned = "alignment/fastq/aligned-{sample}.fastq"
	output:
		"fastq/trimmed/{sample}.trimmed.fastq.gz",
		"alignment/fastq/unaligned-{sample}.fastq.gz",
		"alignment/fastq/aligned-{sample}.fastq.gz"
	shell:"""
	pigz -fk {input.trim}
	pigz -fk {input.aligned}
	pigz -fk {input.unaligned}
	"""
# Do cross-correlation analysis to determine fragment size and write un-normalized wig files for viewing in IGV (Using spp package)
# Install r libraries - spp, Rsamtools and fastcluster prior to running this job
rule cross_correlation:
	input:
		bam = "alignment/{sample}.bam"
	output:
		length = "cross_correlation/fragment_length-{sample}.txt",
		pdf = "cross_correlation/cross_correlation-{sample}.pdf",
		wig_RPM = "wigfiles/{sample}_normRPM.wig",
		wig_SI = "wigfiles/{sample}_normSI.wig",
		read_counts = "reads_stats/{sample}.txt"
	script:
		"scripts/cross-correlation.R"
#Separate S.cerevisiae and S.pombe wig files
rule extract_ScSp_wigfiles:
	input:
		rpm = "wigfiles/{sample}_normRPM.wig",
		si = "wigfiles/{sample}_normSI.wig"
	output:
		Sc_rpm = "wigfiles/Scer/{sample}_Scer_normRPM.wig",
		Sc_si = "wigfiles/Scer/{sample}_Scer_normSI.wig",
		Sp_si = "wigfiles/Spom/{sample}_Spom_normSI.wig"
	shell:"""
	grep "Scer_" {input.rpm} | sed 's/Scer_//g' | sed r <(head -1 {input.rpm}) - > {output.Sc_rpm}
	grep "Scer_" {input.si} | sed 's/Scer_//g' | sed r <(head -1 {input.si}) - > {output.Sc_si}
	grep "Spom_" {input.si} | sed 's/Spom_//g' | sed r <(head -1 {input.si}) - > {output.Sp_si}
	"""
#Combine the Sc/Sp read counts for all samples and generate plots of percentage spike in for all samples
rule plot_spikein_percentage:
	input:
		files = expand("reads_stats/{sample}.txt", sample=SAMPLES)
	output:
		all_txt = "reads_stats/spikein_proportion_all.txt",
		ip_plots = expand("reads_stats/plots/spike_in_proportion_{ip}.{type}", ip=uniq(IP), type = config["imagetype"]["plot_spikein_percentage"])
	params:
		ip = IP,
		samples = SAMPLES,
		plot_type = config["imagetype"]["plot_spikein_percentage"]
	script: "scripts/plot_spikein_proportion.R"




