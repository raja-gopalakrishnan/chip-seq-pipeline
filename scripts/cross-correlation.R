library(spp)
library(ggplot2)
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
p = {
	par(mar = c(3.5,3.5,1.0,0.5), mgp = c(2,0.65,0), cex = 0.8)
	plot(binding.characteristics$cross.correlation,type='l',xlab="Strand shift",ylab="Cross-correlation", main="");
	abline(v=binding.characteristics$peak$x,lty=2,col=2)
}
ggsave(filename=snakemake@output[["png"]], plot = p,width=5, height=5, units="in")
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
