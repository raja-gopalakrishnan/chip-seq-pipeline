#!/usr/bin/env python

configfile: "config.yaml"
SAMPLES = config["samples"]
ANNOTATIONS = config["annotations"]

PEAK_TREATMENT = list(config["callpeaks"].keys())
PEAK_CONTROL = list(config["callpeaks"].values())
NORMALIZATIONS = list(config["normalizations"].keys())

IP = [v["ip"] for k,v in SAMPLES.items()]
GTYPE = [v["gtype"] for k,v in SAMPLES.items()]
REPLICATE = [v["replicate"] for k,v in SAMPLES.items()]

#Function that outputs unique items in a list in the same order as present in the list
def uniq(input):
    output=[]
    for x in input:
        if x not in output:
            output.append(x)
    return output

#Function that returns the names of samples associated with a particular attribute
#Example: attribute = ip and value = ha will give the names of all samples have been subjected to a ha ip
def names(dictionary, attribute, value):
	output=[]
    #example: If dictionary is SAMPLES, level_one is a list of all sample names
	#example: If dictionary is SAMPLES, level_two is a list of the dictionaries associated with each sample name
	level_one=[k for k,v in dictionary.items()]
	level_two=[v for k,v in dictionary.items()]
	for i in range(0,len(level_two)):
		if level_two[i][attribute]==value :
			output.append(level_one[i])
	return output

localrules: all,
   	    make_barcode_file

rule all:
	input:
		#make_barcode_file
		"barcodes.tsv",
		#fastqc_raw
		"qual_ctrl/raw",
		#fastqc_processed
		expand("qual_ctrl/trimmed-{sample}", sample=SAMPLES),
		#bowtie
		expand("alignment/fastq/unaligned-{sample}.fastq.gz", sample=SAMPLES),
		expand("alignment/fastq/aligned-{sample}.fastq.gz", sample=SAMPLES),
		#gzip_loose_fastq
		expand("fastq/trimmed/{sample}.trimmed.fastq.gz", sample=SAMPLES),
		#samtools_index
		expand("alignment/{sample}{species}.bam.bai", sample = SAMPLES, species = ["","_Scer","_Spom"]),
		#bam_separate_species
		expand("alignment/{sample}_{species}.bam", sample = SAMPLES, species = ["Scer", "Spom"]),
		#cross_correlation
		expand("cross_correlation/fragment_length-{ip}-{gtype}-{replicate}.txt", ip=uniq(IP), gtype=uniq(GTYPE), replicate=uniq(REPLICATE)),
		#bowtie_summary
		"alignment/bowtie_read_number_summary.txt",
		#generate_coverage
		expand("wigfiles/{species}/{sample}_{species}_unnormalized.wig", sample = SAMPLES, species = ["Scer","Spom"]),
		#wig_to_bedgraph
		#expand("bedgraphs/{sample}_{species}.bdg", sample = SAMPLES, species = ["Scer","Spom"]),
		#normalize
		expand("bedgraphs/{species}_normSI/{sample}_{species}_normSI.bdg", sample = SAMPLES, species = ["Scer","Spom"]),
		#plot_spikein_percentage
		expand("reads_stats/{readtype}_spikein_proportion_all.txt", readtype = ["actual","apparent"]),
		expand("plots/{readtype}_spike_in_proportion/{readtype}_spike_in_proportion_{ip}.{type}",readtype = ["actual","apparent"],ip=uniq(IP), type = config["imagetype"]["plot_spikein_percentage"]),
		#make_windows
		expand("bedfiles/{species}_100bp_windows.bed", species = ["Scer","Spom"]),
		#get_window_coverage
		expand("bedgraphs/100bp_windows/{species}_norm{norm}/{sample}_100bp_windows_{species}_norm{norm}.bdg", sample = SAMPLES, species = ["Scer","Spom"], norm = ["RPM", "SI"]),
		#merge_bedgraphs
		expand("bedgraphs/merge_bedgraphs/all_100bp_windows_{species}_norm{norm}.bdg", species = ["Scer","Spom"], norm=["RPM","SI"]),
		#plot_correlations
		expand("plots/correlations/{species}_norm{norm}/{ip}_correlations_{species}_norm{norm}.{type}", ip=uniq(IP), type=config["imagetype"]["plot_correlations"], species=["Scer","Spom"], norm=["RPM", "SI"]),
		#bedgraph_to_bigwig
		expand("bigwigfiles/{species}_norm{norm}/{sample}_{species}_norm{norm}.bw", sample = SAMPLES, norm=["RPM","SI"], species = ["Scer","Spom"]),
		#deeptools_matrix_individual
		expand("matrix/{annotation}/individual/{species}_norm{norm}/{sample}_{species}_norm{norm}_{annotation}.mat.gz", norm = ["RPM", "SI"], sample = SAMPLES, annotation = names(ANNOTATIONS,"species","all"), species = ["Scer","Spom"]),
		expand("matrix/{annotation}/individual/{species}_norm{norm}/{sample}_{species}_norm{norm}_{annotation}.mat.gz", norm = ["RPM", "SI"], sample = SAMPLES, annotation = names(ANNOTATIONS,"species","Scer"), species = ["Scer"]),
		#plot_heatmap_individual
		expand("plots/heatmaps/individual/{annotation}_{species}_norm{norm}/{sample}_{species}_norm{norm}_{annotation}_heatmap." + config["imagetype"]["plot_heatmap"], sample = SAMPLES, annotation = names(ANNOTATIONS,"species","all"), norm = ["RPM","SI"], species = ["Scer","Spom"]),
		expand("plots/heatmaps/individual/{annotation}_{species}_norm{norm}/{sample}_{species}_norm{norm}_{annotation}_heatmap." + config["imagetype"]["plot_heatmap"], sample = SAMPLES, annotation = names(ANNOTATIONS,"species","Scer"), norm = ["RPM","SI"], species = ["Scer"]),
		#deeptools_matrix_group
		expand("matrix/{annotation}/group_ip/{species}_norm{norm}/{ip}_{species}_norm{norm}_{annotation}.mat.gz", norm = ["RPM", "SI"], ip = uniq(IP), annotation = names(ANNOTATIONS,"species","all"), species = ["Scer","Spom"]),
		expand("matrix/{annotation}/group_ip/{species}_norm{norm}/{ip}_{species}_norm{norm}_{annotation}.mat.gz", norm = ["RPM", "SI"], ip = uniq(IP), annotation = names(ANNOTATIONS,"species","Scer"), species = ["Scer"]),
		#plot_heatmap_group
		expand("plots/heatmaps/group_ip/{annotation}_{species}_norm{norm}/{ip}_{species}_norm{norm}_{annotation}_heatmap." + config["imagetype"]["plot_heatmap"], norm = ["RPM", "SI"], ip = uniq(IP), annotation = names(ANNOTATIONS,"species","all"), species = ["Scer","Spom"]),
		expand("plots/heatmaps/group_ip/{annotation}_{species}_norm{norm}/{ip}_{species}_norm{norm}_{annotation}_heatmap." + config["imagetype"]["plot_heatmap"], norm = ["RPM", "SI"], ip = uniq(IP), annotation = names(ANNOTATIONS,"species","Scer"), species = ["Scer"]),
		#macs2
		#expand("peakcalling/{tsample}-{gtype}_{csample}-{gtype}_peaks.narrowPeak", tsample = PEAK_TREATMENT, csample = uniq(PEAK_CONTROL), gtype = uniq(GTYPE)),
		#separate_peaks
		expand("peakcalling/{species}/{tsample}-{gtype}_{csample}-{gtype}_peaks.narrowPeak", species = ["Scer","Spom"], tsample = PEAK_TREATMENT, csample = uniq(PEAK_CONTROL), gtype = uniq(GTYPE)),
		#gzip_deeptools_matrix_individual
		expand("matrix/{annotation}/individual/{species}_norm{norm}/{sample}_{species}_norm{norm}_{annotation}.tsv.gz", norm = ["RPM", "SI"], sample = SAMPLES, annotation = names(ANNOTATIONS,"species","all"), species = ["Scer","Spom"]),
		expand("matrix/{annotation}/individual/{species}_norm{norm}/{sample}_{species}_norm{norm}_{annotation}.tsv.gz", norm = ["RPM", "SI"], sample = SAMPLES, annotation = names(ANNOTATIONS,"species","Scer"), species = ["Scer"]),
		#gzip_deeptools_matrix_group
		expand("matrix/{annotation}/group_ip/{species}_norm{norm}/{ip}_{species}_norm{norm}_{annotation}.tsv.gz", norm = ["RPM", "SI"], ip = uniq(IP), annotation = names(ANNOTATIONS,"species","all"), species = ["Scer","Spom"]),
		expand("matrix/{annotation}/group_ip/{species}_norm{norm}/{ip}_{species}_norm{norm}_{annotation}.tsv.gz", norm = ["RPM", "SI"], ip = uniq(IP), annotation = names(ANNOTATIONS,"species","Scer"), species = ["Scer"]),
		#multi_metagene
		expand("matrix/group_ip_melted/{annotation}/{species}_norm{norm}/{ip}_{species}_norm{norm}_{annotation}_meltedmat.txt", norm = ["RPM", "SI"], ip = uniq(IP), annotation = names(ANNOTATIONS,"species","all"), species = ["Scer","Spom"]),
		expand("matrix/group_ip_melted/{annotation}/{species}_norm{norm}/{ip}_{species}_norm{norm}_{annotation}_meltedmat.txt", norm = ["RPM", "SI"], ip = uniq(IP), annotation = names(ANNOTATIONS,"species","Scer"), species = ["Scer"]),
		#multi_metagene_divbylib
		expand("dummyfile/multi_metagene_divbylib_{ip}_{species}_norm{norm}_{annotation}.txt", norm = ["RPM", "SI"], ip = NORMALIZATIONS, annotation = names(ANNOTATIONS,"species","all"), species = ["Scer","Spom"]),
		expand("dummyfile/multi_metagene_divbylib_{ip}_{species}_norm{norm}_{annotation}.txt", norm = ["RPM", "SI"], ip = NORMALIZATIONS, annotation = names(ANNOTATIONS,"species","Scer"), species = ["Scer"])

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
# Make separate BAM files for S.cerevisiae and S.pombe
rule bam_separate_species:
    input:
        "alignment/{sample}.bam",
        "alignment/{sample}.bam.bai"
    output:
        "alignment/{sample}_{species}.bam"
    params:
    	"genome/combined_genome.chrom.sizes"
    log: "logs/bam_separate_species/bam_separate_species-{sample}_{species}.log"
    shell: """
        (samtools view -b {input} $(grep {wildcards.species} {params} | awk 'BEGIN{{FS="\t"; ORS=" "}}{{print $1}}') > {output}) &> {log}
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
# Do cross-correlation analysis to determine fragment size and write normalized wig files for viewing in IGV (Using spp package)
# Install r libraries - spp, Rsamtools and fastcluster prior to running this job
rule cross_correlation:
	input:
		chip_bam = "alignment/{ip}-{gtype}-{replicate}.bam",
		control_bam = "alignment/input-{gtype}-{replicate}.bam"
	output:
		length = "cross_correlation/fragment_length-{ip}-{gtype}-{replicate}.txt",
		png = "plots/cross_correlation/cross_correlation-{ip}-{gtype}-{replicate}.png",
		read_counts = "reads_stats/{ip}-{gtype}-{replicate}_actualreadcounts.txt",
		apparent_read_counts = "reads_stats/{ip}-{gtype}-{replicate}_apparentreadcounts.txt"
	wildcard_constraints:
		ip = "[a-z0-9]+",
		replicate = "\d+"
	script:
		"scripts/cross-correlation.R"
#Make coverage wigfile
rule generate_coverage:
	input:
		bam = "alignment/{sample}_{species}.bam",
		chrsizes = "genome/{species}.chrom.sizes",
		fragment = "cross_correlation/fragment_length-{sample}.txt",
		index = "alignment/{sample}_{species}.bam.bai"
	output:
		wig = "wigfiles/{species}/{sample}_{species}_unnormalized.wig"
	log: "logs/generate_coverage/generate_coverage_{sample}_{species}.log"
	shell:"""
	read_length=$(samtools view {input.bam} | awk 'BEGIN{{sum=0}}{{sum+=length($10)}}END{{printf "%.0f",sum/FNR}}')
	(igvtools count -w 20 -e $(($(less {input.fragment})-$read_length)) {input.bam} {output.wig} {input.chrsizes}) &> {log}
	"""
#Convert wig to bedgraph
rule wig_to_bedgraph:
	input:
		"wigfiles/{species}/{sample}_{species}_unnormalized.wig"
	output:
		temp("bedgraphs/{sample}_{species}.bdg")
	shell:"""
	wig2bed <{input} | awk 'BEGIN{{FS=OFS="\t"}}{{print $1,$2,$3,$5}}' | sed 's/{wildcards.species}_//g' | sort -k1,1 -k2,2n - > {output}
	"""
#Normalize by spike-in or RPM
rule normalize:
	input:
		bedgraph = "bedgraphs/{sample}_{species}.bdg",
		normfactor = "reads_stats/{sample}_apparentreadcounts.txt"
	output:
		rpm = "bedgraphs/{species}_normRPM/{sample}_{species}_normRPM.bdg",
		si = "bedgraphs/{species}_normSI/{sample}_{species}_normSI.bdg"
	shell:"""
	 awk -v norm=$(sed -n '1p' {input.normfactor}) 'BEGIN{{FS=OFS="\t"}}{{print $1,$2,$3,$4/(norm/1000000)}}' {input.bedgraph} > {output.rpm}
	 awk -v norm=$(sed -n '2p' {input.normfactor}) 'BEGIN{{FS=OFS="\t"}}{{print $1,$2,$3,$4/(norm/1000000)}}' {input.bedgraph} > {output.si}
	"""
#Combine the Sc/Sp read counts for all samples and generate plots of percentage spike in for all samples
#Install R libraries ggplot2, reshape2, Cairo and dplyr prior to running this job
rule plot_spikein_percentage:
	input:
		files = expand("reads_stats/{sample}_{{readtype}}readcounts.txt", sample=SAMPLES)
	output:
		all_txt = "reads_stats/{readtype}_spikein_proportion_all.txt",
		ip_plots = expand("plots/{{readtype}}_spike_in_proportion/{{readtype}}_spike_in_proportion_{ip}.{type}", ip=uniq(IP), type = config["imagetype"]["plot_spikein_percentage"])
	params:
		ip = IP,
		samples = SAMPLES,
		plot_type = config["imagetype"]["plot_spikein_percentage"]
	script: "scripts/plot_spikein_proportion.R"
#Make a bed file with 100 bp windows for cerevisiae, pombe and combined genomes
rule make_windows:
	input:
		"genome/{species}.tsv"
	output:
		"bedfiles/{species}_100bp_windows.bed"
	shell:"""
	bedtools makewindows -g {input} -w 100 | sort -k1,1 -k2,2n - > {output}
	"""
#Get a bedgraph of coverage by mapping to a regularly binned genome file. This can then by used for correlations and input normalizations.
rule get_window_coverage:
	input:
		combined_bed = "bedfiles/{species}_100bp_windows.bed",
		bedgraph = "bedgraphs/{species}_norm{norm}/{sample}_{species}_norm{norm}.bdg"
	output:
		"bedgraphs/100bp_windows/{species}_norm{norm}/{sample}_100bp_windows_{species}_norm{norm}.bdg"
	shell:"""
	bedtools map -a {input.combined_bed} -b {input.bedgraph} -c 4 -o mean -null 0 | sort -k1,1 -k2,2n - > {output}
	"""
#Merge all the bedgraphs of all samples making it easy to plot correlations
rule merge_bedgraphs:
	input:
		files = expand("bedgraphs/100bp_windows/{{species}}_norm{{norm}}/{sample}_100bp_windows_{{species}}_norm{{norm}}.bdg", sample = SAMPLES)		
	output:
		"bedgraphs/merge_bedgraphs/all_100bp_windows_{species}_norm{norm}.bdg"
	params:
		names = expand("{sample}", sample = SAMPLES)
	shell:"""
	bedtools unionbedg -i {input.files} -header -names {params.names} > {output}
	"""
#Plot correlation profiles for spike-in normalized datasets
rule plot_correlations:
	input:
		"bedgraphs/merge_bedgraphs/all_100bp_windows_{species}_norm{norm}.bdg"
	output:
		"plots/correlations/{species}_norm{norm}/{ip}_correlations_{species}_norm{norm}.{type}"
	params:
		samplelist = lambda wildcards: names(SAMPLES,"ip",wildcards.ip),
		pcount = 0.1
	script: "scripts/plotcorr.R"
#Convert bedgraph to bigwig
rule bedgraph_to_bigwig:
	input:
		bedgraph = "bedgraphs/{species}_norm{norm}/{sample}_{species}_norm{norm}.bdg",
		chrsizes = "genome/{species}.tsv"
	output:
		"bigwigfiles/{species}_norm{norm}/{sample}_{species}_norm{norm}.bw"
	shell:"""
	bedGraphToBigWig {input.bedgraph} {input.chrsizes} {output}
	"""
# Generate matrices for plotting heatmaps
rule deeptools_matrix_individual:
	input:
		annotation = lambda wildcards: "genome/annotations/"+ wildcards.species +"_" + config["annotations"][wildcards.annotation]["bed_suffix"] + ".bed",
		bigwig = "bigwigfiles/{species}_norm{norm}/{sample}_{species}_norm{norm}.bw"
	output:
		dtfile = "matrix/{annotation}/individual/{species}_norm{norm}/{sample}_{species}_norm{norm}_{annotation}.mat.gz",
		matrix = temp("matrix/{annotation}/individual/{species}_norm{norm}/{sample}_{species}_norm{norm}_{annotation}.tsv")
	params:
		refpoint = lambda wildcards: config["annotations"][wildcards.annotation]["refpoint"],
		upstream = lambda wildcards: config["annotations"][wildcards.annotation]["upstream"],
		downstream = lambda wildcards: config["annotations"][wildcards.annotation]["downstream"],
		binsize = lambda wildcards: config["annotations"][wildcards.annotation]["binsize"],
		sort = lambda wildcards: config["annotations"][wildcards.annotation]["sort"],
		sortusing = lambda wildcards: config["annotations"][wildcards.annotation]["sortby"],
		binstat = lambda wildcards: config["annotations"][wildcards.annotation]["binstat"],
		bodylength = lambda wildcards: config["annotations"][wildcards.annotation]["bodylength"],
		nan = lambda wildcards: config["annotations"][wildcards.annotation]["nan_afterend"]
	threads: config["threads"]
	log: "logs/deeptools_matrix/{sample}_norm{norm}_{annotation}.log"
	run:
		if config["annotations"][wildcards.annotation]["reference_point"]=="y":
			shell("(computeMatrix reference-point -R {input.annotation} -S {input.bigwig} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.downstream} --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads} {params.nan}) &> {log}")
		else:
			shell("(computeMatrix scale-regions -R {input.annotation} -S {input.bigwig} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.downstream} --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads} --regionBodyLength {params.bodylength}) &> {log}")
# Generate a heatmap using the deeptools matrix
rule plot_heatmap_individual:
	input: 
		dtfile = "matrix/{annotation}/individual/{species}_norm{norm}/{sample}_{species}_norm{norm}_{annotation}.mat.gz"
	output:
		plot = "plots/heatmaps/individual/{annotation}_{species}_norm{norm}/{sample}_{species}_norm{norm}_{annotation}_heatmap." + config["imagetype"]["plot_heatmap"]
	params:
		sort = lambda wildcards: config["annotations"][wildcards.annotation]["sort"],
		sortusing = lambda wildcards: config["annotations"][wildcards.annotation]["sortby"],
		missingdatacolor = lambda wildcards: config["annotations"][wildcards.annotation]["missingdatacolor"]
	shell:"""
	plotHeatmap -m {input.dtfile} -out {output.plot} --sortRegions {params.sort} --sortUsing {params.sortusing} --missingDataColor {params.missingdatacolor}
	"""
# Generate matrices for plotting heatmaps by IP
rule deeptools_matrix_group:
	input:
		annotation = lambda wildcards: "genome/annotations/"+ wildcards.species +"_" + config["annotations"][wildcards.annotation]["bed_suffix"] + ".bed",
		bigwig = lambda wildcards: expand("bigwigfiles/{species}_norm{norm}/{sample}_{species}_norm{norm}.bw", sample = names(SAMPLES,"ip",wildcards.ip), species = wildcards.species, norm=wildcards.norm)
	output:
		dtfile = "matrix/{annotation}/group_ip/{species}_norm{norm}/{ip}_{species}_norm{norm}_{annotation}.mat.gz",
		matrix = temp("matrix/{annotation}/group_ip/{species}_norm{norm}/{ip}_{species}_norm{norm}_{annotation}.tsv")
	params:
		refpoint = lambda wildcards: config["annotations"][wildcards.annotation]["refpoint"],
		upstream = lambda wildcards: config["annotations"][wildcards.annotation]["upstream"],
		downstream = lambda wildcards: config["annotations"][wildcards.annotation]["downstream"],
		binsize = lambda wildcards: config["annotations"][wildcards.annotation]["binsize"],
		sort = lambda wildcards: config["annotations"][wildcards.annotation]["sort"],
		sortusing = lambda wildcards: config["annotations"][wildcards.annotation]["sortby"],
		binstat = lambda wildcards: config["annotations"][wildcards.annotation]["binstat"],
		bodylength = lambda wildcards: config["annotations"][wildcards.annotation]["bodylength"],
		nan = lambda wildcards: config["annotations"][wildcards.annotation]["nan_afterend"]
	threads: config["threads"]
	log: "logs/deeptools_matrix/{ip}_norm{norm}_{annotation}.log"
	run:
		if config["annotations"][wildcards.annotation]["reference_point"]=="y":
			shell("(computeMatrix reference-point -R {input.annotation} -S {input.bigwig} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.downstream} --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads} {params.nan}) &> {log}")
		else:
			shell("(computeMatrix scale-regions -R {input.annotation} -S {input.bigwig} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.downstream} --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads} --regionBodyLength {params.bodylength}) &> {log}")
#Generate heatmap using plot_heatmap
rule plot_heatmap_group:
	input:
		dtfile = "matrix/{annotation}/group_ip/{species}_norm{norm}/{ip}_{species}_norm{norm}_{annotation}.mat.gz"
	output:
		plot = "plots/heatmaps/group_ip/{annotation}_{species}_norm{norm}/{ip}_{species}_norm{norm}_{annotation}_heatmap." + config["imagetype"]["plot_heatmap"]
	params:
		sort = lambda wildcards: config["annotations"][wildcards.annotation]["sort"],
		sortusing = lambda wildcards: config["annotations"][wildcards.annotation]["sortby"],
		missingdatacolor = lambda wildcards: config["annotations"][wildcards.annotation]["missingdatacolor"]
	shell:"""
	plotHeatmap -m {input.dtfile} -out {output.plot} --sortRegions {params.sort} --sortUsing {params.sortusing} --missingDataColor {params.missingdatacolor}
	"""
#Call peaks
#Makes the assumption that fragment sizes for replicate samples are approximately equal, so extends reads by fragment size of replicate 1
rule macs2:
	input:
		treatment = lambda wildcards: expand("alignment/{sample}.bam", sample = names(SAMPLES,"group",wildcards.tgroup)),
		control = lambda wildcards: expand("alignment/{sample}.bam", sample = names(SAMPLES,"group",wildcards.cgroup)),
		chrsizes = lambda wildcards: "genome/combined_genome.chrom.sizes",
		fragsize = "cross_correlation/fragment_length-{tgroup}-1.txt"
	output:
		xls = temp("peakcalling/{tgroup}_{cgroup}_peaks.xls"),
		peaks = temp("peakcalling/{tgroup}_{cgroup}_peaks.narrowPeak"),
		summits = temp("peakcalling/{tgroup}_{cgroup}_summits.bed"),
		treat_bg = temp("peakcalling/{tgroup}_{cgroup}_treat_pileup.bdg"),
		cntrl_bg = temp("peakcalling/{tgroup}_{cgroup}_control_lambda.bdg")
	params:
		qscore = config["macs2"]["fdr"],
		prefix = lambda wildcards: wildcards.tgroup + "_" + wildcards.cgroup,
		env = config["macs2"]["env_name"]
	conda:
		"envs/python2.yaml"
	log: "logs/macs2/macs2_{tgroup}_{cgroup}.log"
	shell:"""
	(macs2 callpeak -t {input.treatment} -f BAM -g $(awk '{{sum += $2}} END {{print sum}}' {input.chrsizes}) --bdg --outdir peakcalling/ -n {wildcards.tgroup}_{wildcards.cgroup} -q {params.qscore} -c {input.control} --nomodel --call-summits --SPMR --extsize $(less {input.fragsize})) &> {log}
    """
#Separate Scer and Spom peaks
rule separate_peaks:
	input:
		xls = "peakcalling/{tgroup}_{cgroup}_peaks.xls",
		peaks = "peakcalling/{tgroup}_{cgroup}_peaks.narrowPeak",
		summits = "peakcalling/{tgroup}_{cgroup}_summits.bed",
		treat_bg = "peakcalling/{tgroup}_{cgroup}_treat_pileup.bdg",
		cntrl_bg = "peakcalling/{tgroup}_{cgroup}_control_lambda.bdg"
	output:
		xls = "peakcalling/{species}/{tgroup}_{cgroup}_peaks.xls",
		peaks = "peakcalling/{species}/{tgroup}_{cgroup}_peaks.narrowPeak",
		summits = "peakcalling/{species}/{tgroup}_{cgroup}_summits.bed",
		treat_bg = "peakcalling/{species}/{tgroup}_{cgroup}_treat_pileup.bdg",
		cntrl_bg = "peakcalling/{species}/{tgroup}_{cgroup}_control_lambda.bdg"
	shell:"""
	set +e
	grep {wildcards.species} {input.xls} | sed -e 's/{wildcards.tgroup}_{wildcards.cgroup}/{wildcards.tgroup}/g' -e 's/{wildcards.species}_//g' > {output.xls}
	grep {wildcards.species} {input.peaks} | sed -e 's/{wildcards.tgroup}_{wildcards.cgroup}/{wildcards.tgroup}/g' -e 's/{wildcards.species}_//g' > {output.peaks}
	grep {wildcards.species} {input.summits} | sed -e 's/{wildcards.tgroup}_{wildcards.cgroup}/{wildcards.tgroup}/g' -e 's/{wildcards.species}_//g' > {output.summits}
	grep {wildcards.species} {input.treat_bg} | sed 's/{wildcards.species}_//g' > {output.treat_bg}
	grep {wildcards.species} {input.cntrl_bg} | sed 's/{wildcards.species}_//g' > {output.cntrl_bg}
	"""
#Gzip the matrices generated by deeptools
rule gzip_deeptools_matrix_individual:
	input:
		"matrix/{annotation}/individual/{species}_norm{norm}/{sample}_{species}_norm{norm}_{annotation}.tsv"
	output:
		"matrix/{annotation}/individual/{species}_norm{norm}/{sample}_{species}_norm{norm}_{annotation}.tsv.gz"
	shell:"""
	pigz -f {input}
	"""
rule gzip_deeptools_matrix_group:
	input:
		"matrix/{annotation}/group_ip/{species}_norm{norm}/{ip}_{species}_norm{norm}_{annotation}.tsv"
	output:
		"matrix/{annotation}/group_ip/{species}_norm{norm}/{ip}_{species}_norm{norm}_{annotation}.tsv.gz"
	shell:"""
	pigz -f {input}
	"""
#Plot metagenes for all samples subjected to the same IP on the same plot
rule multi_metagene:
	input:
		mat = "matrix/{annotation}/group_ip/{species}_norm{norm}/{ip}_{species}_norm{norm}_{annotation}.tsv.gz"
	output:
		matname = "matrix/group_ip_melted/{annotation}/{species}_norm{norm}/{ip}_{species}_norm{norm}_{annotation}_meltedmat.txt",
		plotname = "plots/metagene/group_ip_div_none/{annotation}/{species}_norm{norm}/{ip}_{species}_norm{norm}_{annotation}_metagene." + config["imagetype"]["multi_metagene"]
	params:
		binsize = lambda wildcards: config["annotations"][wildcards.annotation]["binsize"],
		upstream = lambda wildcards: config["annotations"][wildcards.annotation]["upstream"],
		downstream = lambda wildcards: config["annotations"][wildcards.annotation]["downstream"],
		bodylength = lambda wildcards: config["annotations"][wildcards.annotation]["bodylength"],
		align = lambda wildcards: config["annotations"][wildcards.annotation]["refpoint"]
	script: "scripts/multi_metagene.R"
#Plot metagenes for all samples normalizing to input or another IP
rule multi_metagene_divbylib:
	input:
		expand("matrix/group_ip_melted/{{annotation}}/{{species}}_norm{{norm}}/{ip}_{{species}}_norm{{norm}}_{{annotation}}_meltedmat.txt", ip = uniq(IP)),
		mat_ip = "matrix/group_ip_melted/{annotation}/{species}_norm{norm}/{ip}_{species}_norm{norm}_{annotation}_meltedmat.txt"
	output:
		"dummyfile/multi_metagene_divbylib_{ip}_{species}_norm{norm}_{annotation}.txt"
	params:
		prefix = "matrix/group_ip_melted/{annotation}/{species}_norm{norm}/",
		suffix = "_{species}_norm{norm}_{annotation}_meltedmat.txt",
		control = lambda wildcards: config["normalizations"][wildcards.ip],
		binsize = lambda wildcards: config["annotations"][wildcards.annotation]["binsize"],
		upstream = lambda wildcards: config["annotations"][wildcards.annotation]["upstream"],
		downstream = lambda wildcards: config["annotations"][wildcards.annotation]["downstream"],
		bodylength = lambda wildcards: config["annotations"][wildcards.annotation]["bodylength"],
		align = lambda wildcards: config["annotations"][wildcards.annotation]["refpoint"],
		exclude = config["exclude"]["multi_metagene_divbylib"],
		outpath = "plots/metagene/group_ip_div_lib/{annotation}/{species}_norm{norm}/{ip}_{species}_norm{norm}_{annotation}_div_",
		imagetype = config["imagetype"]["multi_metagene_divbylib"],
		outdir = "plots/metagene/group_ip_div_lib/{annotation}/{species}_norm{norm}/"
	script: "scripts/metagene_divbylib.R"