
######## Snakemake header ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character"
    )
)
snakemake <- Snakemake(
    input = list('alignment/h3k36me3-spt6-1.bam', 'alignment/input-spt6-1.bam', "chip_bam" = 'alignment/h3k36me3-spt6-1.bam', "control_bam" = 'alignment/input-spt6-1.bam'),
    output = list('cross_correlation/fragment_length-h3k36me3-spt6-1.txt', 'plots/cross_correlation/cross_correlation-h3k36me3-spt6-1.png', 'reads_stats/h3k36me3-spt6-1_actualreadcounts.txt', 'reads_stats/h3k36me3-spt6-1_apparentreadcounts.txt', "length" = 'cross_correlation/fragment_length-h3k36me3-spt6-1.txt', "png" = 'plots/cross_correlation/cross_correlation-h3k36me3-spt6-1.png', "read_counts" = 'reads_stats/h3k36me3-spt6-1_actualreadcounts.txt', "apparent_read_counts" = 'reads_stats/h3k36me3-spt6-1_apparentreadcounts.txt'),
    params = list(),
    wildcards = list('h3k36me3', 'spt6', '1', "ip" = 'h3k36me3', "gtype" = 'spt6', "replicate" = '1'),
    threads = 1,
    log = list(),
    resources = list(),
    config = list("samples" = list("ha-wt-1" = list("barcode" = 'ATCACG', "group" = 'ha-wt', "ip" = 'ha', "gtype" = 'wt', "replicate" = 1), "ha-spt6-1" = list("barcode" = 'CGATGT', "group" = 'ha-spt6', "ip" = 'ha', "gtype" = 'spt6', "replicate" = 1), "ha-spt6-set2-1" = list("barcode" = 'TTAGGC', "group" = 'ha-spt6-set2', "ip" = 'ha', "gtype" = 'spt6-set2', "replicate" = 1), "ha-set2-1" = list("barcode" = 'TGACCA', "group" = 'ha-set2', "ip" = 'ha', "gtype" = 'set2', "replicate" = 1), "ha-wt-2" = list("barcode" = 'ACAGTG', "group" = 'ha-wt', "ip" = 'ha', "gtype" = 'wt', "replicate" = 2), "ha-spt6-2" = list("barcode" = 'GCCAAT', "group" = 'ha-spt6', "ip" = 'ha', "gtype" = 'spt6', "replicate" = 2), "ha-spt6-set2-2" = list("barcode" = 'CAGATC', "group" = 'ha-spt6-set2', "ip" = 'ha', "gtype" = 'spt6-set2', "replicate" = 2), "ha-set2-2" = list("barcode" = 'ACTTGA', "group" = 'ha-set2', "ip" = 'ha', "gtype" = 'set2', "replicate" = 2), "pol2-wt-1" = list("barcode" = 'GATCAG', "group" = 'pol2-wt', "ip" = 'pol2', "gtype" = 'wt', "replicate" = 1), "pol2-spt6-1" = list("barcode" = 'TAGCTT', "group" = 'pol2-spt6', "ip" = 'pol2', "gtype" = 'spt6', "replicate" = 1), "pol2-spt6-set2-1" = list("barcode" = 'GGCTAC', "group" = 'pol2-spt6-set2', "ip" = 'pol2', "gtype" = 'spt6-set2', "replicate" = 1), "pol2-set2-1" = list("barcode" = 'CTTGTA', "group" = 'pol2-set2', "ip" = 'pol2', "gtype" = 'set2', "replicate" = 1), "pol2-wt-2" = list("barcode" = 'ATATAGGA', "group" = 'pol2-wt', "ip" = 'pol2', "gtype" = 'wt', "replicate" = 2), "pol2-spt6-2" = list("barcode" = 'AACCGTGT', "group" = 'pol2-spt6', "ip" = 'pol2', "gtype" = 'spt6', "replicate" = 2), "pol2-spt6-set2-2" = list("barcode" = 'AGGTCAGT', "group" = 'pol2-spt6-set2', "ip" = 'pol2', "gtype" = 'spt6-set2', "replicate" = 2), "pol2-set2-2" = list("barcode" = 'CTCTGTCT', "group" = 'pol2-set2', "ip" = 'pol2', "gtype" = 'set2', "replicate" = 2), "h3k36me2-wt-1" = list("barcode" = 'CCATACAC', "group" = 'h3k36me2-wt', "ip" = 'h3k36me2', "gtype" = 'wt', "replicate" = 1), "h3k36me2-spt6-1" = list("barcode" = 'CGCATTAA', "group" = 'h3k36me2-spt6', "ip" = 'h3k36me2', "gtype" = 'spt6', "replicate" = 1), "h3k36me2-spt6-set2-1" = list("barcode" = 'GTCTACAT', "group" = 'h3k36me2-spt6-set2', "ip" = 'h3k36me2', "gtype" = 'spt6-set2', "replicate" = 1), "h3k36me2-set2-1" = list("barcode" = 'GAGTTAAC', "group" = 'h3k36me2-set2', "ip" = 'h3k36me2', "gtype" = 'set2', "replicate" = 1), "h3k36me2-wt-2" = list("barcode" = 'GCAGCCTC', "group" = 'h3k36me2-wt', "ip" = 'h3k36me2', "gtype" = 'wt', "replicate" = 2), "h3k36me2-spt6-2" = list("barcode" = 'TCGCGTAC', "group" = 'h3k36me2-spt6', "ip" = 'h3k36me2', "gtype" = 'spt6', "replicate" = 2), "h3k36me2-spt6-set2-2" = list("barcode" = 'TATACCGT', "group" = 'h3k36me2-spt6-set2', "ip" = 'h3k36me2', "gtype" = 'spt6-set2', "replicate" = 2), "h3k36me2-set2-2" = list("barcode" = 'TGCGGTTA', "group" = 'h3k36me2-set2', "ip" = 'h3k36me2', "gtype" = 'set2', "replicate" = 2), "h3k36me3-wt-1" = list("barcode" = 'AACACCTAC', "group" = 'h3k36me3-wt', "ip" = 'h3k36me3', "gtype" = 'wt', "replicate" = 1), "h3k36me3-spt6-1" = list("barcode" = 'CCTTTACAG', "group" = 'h3k36me3-spt6', "ip" = 'h3k36me3', "gtype" = 'spt6', "replicate" = 1), "h3k36me3-spt6-set2-1" = list("barcode" = 'GGTCCTTGA', "group" = 'h3k36me3-spt6-set2', "ip" = 'h3k36me3', "gtype" = 'spt6-set2', "replicate" = 1), "h3k36me3-set2-1" = list("barcode" = 'TTGAGTGT', "group" = 'h3k36me3-set2', "ip" = 'h3k36me3', "gtype" = 'set2', "replicate" = 1), "h3k36me3-wt-2" = list("barcode" = 'ACTAACTGC', "group" = 'h3k36me3-wt', "ip" = 'h3k36me3', "gtype" = 'wt', "replicate" = 2), "h3k36me3-spt6-2" = list("barcode" = 'CAGGAGGCG', "group" = 'h3k36me3-spt6', "ip" = 'h3k36me3', "gtype" = 'spt6', "replicate" = 2), "h3k36me3-spt6-set2-2" = list("barcode" = 'GTTGTCCCA', "group" = 'h3k36me3-spt6-set2', "ip" = 'h3k36me3', "gtype" = 'spt6-set2', "replicate" = 2), "h3k36me3-set2-2" = list("barcode" = 'TGACGCAT', "group" = 'h3k36me3-set2', "ip" = 'h3k36me3', "gtype" = 'set2', "replicate" = 2), "h3-wt-1" = list("barcode" = 'ATCGCCAGC', "group" = 'h3-wt', "ip" = 'h3', "gtype" = 'wt', "replicate" = 1), "h3-spt6-1" = list("barcode" = 'CATTCCAAG', "group" = 'h3-spt6', "ip" = 'h3', "gtype" = 'spt6', "replicate" = 1), "h3-spt6-set2-1" = list("barcode" = 'GCAAGTAGA', "group" = 'h3-spt6-set2', "ip" = 'h3', "gtype" = 'spt6-set2', "replicate" = 1), "h3-set2-1" = list("barcode" = 'TGATCCGA', "group" = 'h3-set2', "ip" = 'h3', "gtype" = 'set2', "replicate" = 1), "h3-wt-2" = list("barcode" = 'ACGTAGCTC', "group" = 'h3-wt', "ip" = 'h3', "gtype" = 'wt', "replicate" = 2), "h3-spt6-2" = list("barcode" = 'CGAACTGTG', "group" = 'h3-spt6', "ip" = 'h3', "gtype" = 'spt6', "replicate" = 2), "h3-spt6-set2-2" = list("barcode" = 'TAGCTAGTA', "group" = 'h3-spt6-set2', "ip" = 'h3', "gtype" = 'spt6-set2', "replicate" = 2), "h3-set2-2" = list("barcode" = 'GTGGGATA', "group" = 'h3-set2', "ip" = 'h3', "gtype" = 'set2', "replicate" = 2), "input-wt-1" = list("barcode" = 'ATCCTATTC', "group" = 'input-wt', "ip" = 'input', "gtype" = 'wt', "replicate" = 1), "input-spt6-1" = list("barcode" = 'CGGACGTGG', "group" = 'input-spt6', "ip" = 'input', "gtype" = 'spt6', "replicate" = 1), "input-spt6-set2-1" = list("barcode" = 'GCGTTTCGA', "group" = 'input-spt6-set2', "ip" = 'input', "gtype" = 'spt6-set2', "replicate" = 1), "input-set2-1" = list("barcode" = 'TATCTCCG', "group" = 'input-set2', "ip" = 'input', "gtype" = 'set2', "replicate" = 1), "input-wt-2" = list("barcode" = 'CACAGTTGG', "group" = 'input-wt', "ip" = 'input', "gtype" = 'wt', "replicate" = 2), "input-spt6-2" = list("barcode" = 'GTGACTACA', "group" = 'input-spt6', "ip" = 'input', "gtype" = 'spt6', "replicate" = 2), "input-spt6-set2-2" = list("barcode" = 'TGAGAGTG', "group" = 'input-spt6-set2', "ip" = 'input', "gtype" = 'spt6-set2', "replicate" = 2), "input-set2-2" = list("barcode" = 'AATGCTGAC', "group" = 'input-set2', "ip" = 'input', "gtype" = 'set2', "replicate" = 2)), "fastq" = 'Undetermined_S0.R1.fastq.gz', "demultiplex" = list("mismatch" = 1), "cutadapt" = list("qual_cutoff" = 20, "min_length" = 5), "threads" = 4, "imagetype" = list("plot_spikein_percentage" = 'png', "plot_correlations" = 'png', "plot_heatmap" = 'png'), "normalizations" = list("ha" = 'pol2', "pol2" = 'input', "h3k36me3" = 'h3', "h3k36me2" = 'h3', "h3" = 'input'), "annotations" = list("align-tss" = list("refpoint" = 'TSS', "upstream" = 300, "downstream" = 4000, "binstat" = 'mean', "binsize" = 20, "nan_afterend" = 'y', "sort" = 'ascend', "sortby" = 'region_length', "missingdatacolor" = 1, "reference_point" = 'y', "bodylength" = '', "bed_suffix" = 'polIItranscripts-adjustedTSS'), "align-tes" = list("refpoint" = 'TES', "upstream" = 4000, "downstream" = 300, "binstat" = 'mean', "binsize" = 20, "nan_afterend" = 'y', "sort" = 'ascend', "sortby" = 'region_length', "missingdatacolor" = 1, "reference_point" = 'y', "bodylength" = '', "bed_suffix" = 'polIItranscripts-adjustedTSS'), "align-tss-tes" = list("refpoint" = 'TSS-TES', "upstream" = 500, "downstream" = 500, "binstat" = 'mean', "binsize" = 20, "nan_afterend" = 'y', "sort" = 'ascend', "sortby" = 'region_length', "missingdatacolor" = 1, "reference_point" = 'n', "bodylength" = 2000, "bed_suffix" = 'polIItranscripts-adjustedTSS')), "callpeaks" = list("h3k36me3" = 'input', "h3k36me2" = 'input', "ha" = 'input', "pol2" = 'input', "h3" = 'input'), "macs2" = list("fdr" = 0.05, "env_name" = 'python2')),
    rule = 'cross_correlation'
)
######## Original script #########
library(spp)
library(Cairo)
options(scipen=100)
#Read in bam file
chip.data = read.bam.tags(snakemake@input[["chip_bam"]])
control.data = read.bam.tags(snakemake@input[["control_bam"]])

# Calculate cross-correlation profile
binding.characteristics = get.binding.characteristics(chip.data, srange=c(0,500), bin=5, remove.tag.anomalies=T, accept.all.tags=T)

#Get fragment size
frag_size = binding.characteristics$peak$x
write.table(frag_size, snakemake@output[["length"]], quote=F, row.names=F, col.names=F)
tag.shift = frag_size/2

#Make plot showing cross-correlation profile
Cairo(file=snakemake@output[["png"]],type="png",width=5,height=5,units="in",dpi=300)
par(mar = c(3.5,3.5,1.0,0.5), mgp = c(2,0.65,0), cex = 0.8)
plot(binding.characteristics$cross.correlation,type='l',xlab="Strand shift",ylab="Cross-correlation", main="");
abline(v=binding.characteristics$peak$x,lty=2,col=2)
dev.off()

#Remove positions with abnormally high number of reads as it probably indicates some technical bias (threshold = 10-fold)
#Threshold fold-over background density above which the position is capped to the maximum statistically likely given local tag density = 4
chip.data=select.informative.tags(chip.data)
chip.data=remove.local.tag.anomalies(chip.data)

control.data=select.informative.tags(control.data)
control.data=remove.local.tag.anomalies(control.data)

#Extract S.cerevisiae and S. pombe read counts.
#Stored as a vector: Position 1 represent chip (IP) reads, Position 2 represents control (input) reads
scer_reads = c(0,0)
spom_reads = c(0,0)
total_reads = c(0,0)

#Get chromosome names
if(sum(names(chip.data)!=names(control.data))==0)
chr_names = names(chip.data)

#Go through the chromosome names in each list and sum the number of reads
for(i in chr_names)
{
	chip.positions = eval(parse(text = paste("chip.data$",i,sep="")))
	control.positions = eval(parse(text = paste("control.data$",i,sep="")))

	organism = substring(i,1,5)
	if(organism == "Scer_")
		scer_reads = scer_reads + c(length(chip.positions),length(control.positions))
	if(organism == "Spom_")
		spom_reads = spom_reads + c(length(chip.positions),length(control.positions))

	total_reads = total_reads + c(length(chip.positions),length(control.positions))
}

#Get spike in and RPM normalization factors
rpm_factor = scer_reads[1]/1000000
si_factor = (spom_reads[1]/1000000)*(0.1*total_reads[2])/spom_reads[2]

#Make a table compiling cerevisiae, pombe and total read counts
readcount_table = c(scer_reads[1],spom_reads[1],total_reads[1])

#Make a table of normalization factors
apparent_readcount_table = c(rpm_factor*1000000, round(si_factor*1000000), (rpm_factor*1000000+round(si_factor*1000000)))

#write all tables
write.table(readcount_table,snakemake@output[["read_counts"]], sep="\t", row.names=F, col.names=F, quote=F)
write.table(apparent_readcount_table,snakemake@output[["apparent_read_counts"]], sep="\t", row.names=F, col.names=F, quote=F)
