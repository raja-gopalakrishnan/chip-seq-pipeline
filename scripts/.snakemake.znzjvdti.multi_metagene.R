
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
    input = list('matrix/align-tss/group_ip/Spom_normSI/pol2_Spom_normSI_align-tss.tsv.gz', "mat" = 'matrix/align-tss/group_ip/Spom_normSI/pol2_Spom_normSI_align-tss.tsv.gz'),
    output = list('matrix/group_ip_melted/align-tss/Spom_normSI/pol2_Spom_normSI_align-tss_meltedmat.txt', 'plots/metagene/group_ip_div_none/align-tss/Spom_normSI/pol2_Spom_normSI_align-tss_metagene.png', "matname" = 'matrix/group_ip_melted/align-tss/Spom_normSI/pol2_Spom_normSI_align-tss_meltedmat.txt', "plotname" = 'plots/metagene/group_ip_div_none/align-tss/Spom_normSI/pol2_Spom_normSI_align-tss_metagene.png'),
    params = list(20, 300, 4000, 0, 'TSS', "binsize" = 20, "upstream" = 300, "downstream" = 4000, "bodylength" = 0, "align" = 'TSS'),
    wildcards = list('align-tss', 'Spom', 'SI', 'pol2', "annotation" = 'align-tss', "species" = 'Spom', "norm" = 'SI', "ip" = 'pol2'),
    threads = 1,
    log = list(),
    resources = list(),
    config = list("samples" = list("ha-wt-1" = list("barcode" = 'ATCACG', "group" = 'ha-wt', "ip" = 'ha', "gtype" = 'wt', "replicate" = 1), "ha-spt6-1" = list("barcode" = 'CGATGT', "group" = 'ha-spt6', "ip" = 'ha', "gtype" = 'spt6', "replicate" = 1), "ha-spt6-set2-1" = list("barcode" = 'TTAGGC', "group" = 'ha-spt6-set2', "ip" = 'ha', "gtype" = 'spt6-set2', "replicate" = 1), "ha-set2-1" = list("barcode" = 'TGACCA', "group" = 'ha-set2', "ip" = 'ha', "gtype" = 'set2', "replicate" = 1), "ha-wt-2" = list("barcode" = 'ACAGTG', "group" = 'ha-wt', "ip" = 'ha', "gtype" = 'wt', "replicate" = 2), "ha-spt6-2" = list("barcode" = 'GCCAAT', "group" = 'ha-spt6', "ip" = 'ha', "gtype" = 'spt6', "replicate" = 2), "ha-spt6-set2-2" = list("barcode" = 'CAGATC', "group" = 'ha-spt6-set2', "ip" = 'ha', "gtype" = 'spt6-set2', "replicate" = 2), "ha-set2-2" = list("barcode" = 'ACTTGA', "group" = 'ha-set2', "ip" = 'ha', "gtype" = 'set2', "replicate" = 2), "pol2-wt-1" = list("barcode" = 'GATCAG', "group" = 'pol2-wt', "ip" = 'pol2', "gtype" = 'wt', "replicate" = 1), "pol2-spt6-1" = list("barcode" = 'TAGCTT', "group" = 'pol2-spt6', "ip" = 'pol2', "gtype" = 'spt6', "replicate" = 1), "pol2-spt6-set2-1" = list("barcode" = 'GGCTAC', "group" = 'pol2-spt6-set2', "ip" = 'pol2', "gtype" = 'spt6-set2', "replicate" = 1), "pol2-set2-1" = list("barcode" = 'CTTGTA', "group" = 'pol2-set2', "ip" = 'pol2', "gtype" = 'set2', "replicate" = 1), "pol2-wt-2" = list("barcode" = 'ATATAGGA', "group" = 'pol2-wt', "ip" = 'pol2', "gtype" = 'wt', "replicate" = 2), "pol2-spt6-2" = list("barcode" = 'AACCGTGT', "group" = 'pol2-spt6', "ip" = 'pol2', "gtype" = 'spt6', "replicate" = 2), "pol2-spt6-set2-2" = list("barcode" = 'AGGTCAGT', "group" = 'pol2-spt6-set2', "ip" = 'pol2', "gtype" = 'spt6-set2', "replicate" = 2), "pol2-set2-2" = list("barcode" = 'CTCTGTCT', "group" = 'pol2-set2', "ip" = 'pol2', "gtype" = 'set2', "replicate" = 2), "h3k36me2-wt-1" = list("barcode" = 'CCATACAC', "group" = 'h3k36me2-wt', "ip" = 'h3k36me2', "gtype" = 'wt', "replicate" = 1), "h3k36me2-spt6-1" = list("barcode" = 'CGCATTAA', "group" = 'h3k36me2-spt6', "ip" = 'h3k36me2', "gtype" = 'spt6', "replicate" = 1), "h3k36me2-spt6-set2-1" = list("barcode" = 'GTCTACAT', "group" = 'h3k36me2-spt6-set2', "ip" = 'h3k36me2', "gtype" = 'spt6-set2', "replicate" = 1), "h3k36me2-set2-1" = list("barcode" = 'GAGTTAAC', "group" = 'h3k36me2-set2', "ip" = 'h3k36me2', "gtype" = 'set2', "replicate" = 1), "h3k36me2-wt-2" = list("barcode" = 'GCAGCCTC', "group" = 'h3k36me2-wt', "ip" = 'h3k36me2', "gtype" = 'wt', "replicate" = 2), "h3k36me2-spt6-2" = list("barcode" = 'TCGCGTAC', "group" = 'h3k36me2-spt6', "ip" = 'h3k36me2', "gtype" = 'spt6', "replicate" = 2), "h3k36me2-spt6-set2-2" = list("barcode" = 'TATACCGT', "group" = 'h3k36me2-spt6-set2', "ip" = 'h3k36me2', "gtype" = 'spt6-set2', "replicate" = 2), "h3k36me2-set2-2" = list("barcode" = 'TGCGGTTA', "group" = 'h3k36me2-set2', "ip" = 'h3k36me2', "gtype" = 'set2', "replicate" = 2), "h3k36me3-wt-1" = list("barcode" = 'AACACCTAC', "group" = 'h3k36me3-wt', "ip" = 'h3k36me3', "gtype" = 'wt', "replicate" = 1), "h3k36me3-spt6-1" = list("barcode" = 'CCTTTACAG', "group" = 'h3k36me3-spt6', "ip" = 'h3k36me3', "gtype" = 'spt6', "replicate" = 1), "h3k36me3-spt6-set2-1" = list("barcode" = 'GGTCCTTGA', "group" = 'h3k36me3-spt6-set2', "ip" = 'h3k36me3', "gtype" = 'spt6-set2', "replicate" = 1), "h3k36me3-set2-1" = list("barcode" = 'TTGAGTGT', "group" = 'h3k36me3-set2', "ip" = 'h3k36me3', "gtype" = 'set2', "replicate" = 1), "h3k36me3-wt-2" = list("barcode" = 'ACTAACTGC', "group" = 'h3k36me3-wt', "ip" = 'h3k36me3', "gtype" = 'wt', "replicate" = 2), "h3k36me3-spt6-2" = list("barcode" = 'CAGGAGGCG', "group" = 'h3k36me3-spt6', "ip" = 'h3k36me3', "gtype" = 'spt6', "replicate" = 2), "h3k36me3-spt6-set2-2" = list("barcode" = 'GTTGTCCCA', "group" = 'h3k36me3-spt6-set2', "ip" = 'h3k36me3', "gtype" = 'spt6-set2', "replicate" = 2), "h3k36me3-set2-2" = list("barcode" = 'TGACGCAT', "group" = 'h3k36me3-set2', "ip" = 'h3k36me3', "gtype" = 'set2', "replicate" = 2), "h3-wt-1" = list("barcode" = 'ATCGCCAGC', "group" = 'h3-wt', "ip" = 'h3', "gtype" = 'wt', "replicate" = 1), "h3-spt6-1" = list("barcode" = 'CATTCCAAG', "group" = 'h3-spt6', "ip" = 'h3', "gtype" = 'spt6', "replicate" = 1), "h3-spt6-set2-1" = list("barcode" = 'GCAAGTAGA', "group" = 'h3-spt6-set2', "ip" = 'h3', "gtype" = 'spt6-set2', "replicate" = 1), "h3-set2-1" = list("barcode" = 'TGATCCGA', "group" = 'h3-set2', "ip" = 'h3', "gtype" = 'set2', "replicate" = 1), "h3-wt-2" = list("barcode" = 'ACGTAGCTC', "group" = 'h3-wt', "ip" = 'h3', "gtype" = 'wt', "replicate" = 2), "h3-spt6-2" = list("barcode" = 'CGAACTGTG', "group" = 'h3-spt6', "ip" = 'h3', "gtype" = 'spt6', "replicate" = 2), "h3-spt6-set2-2" = list("barcode" = 'TAGCTAGTA', "group" = 'h3-spt6-set2', "ip" = 'h3', "gtype" = 'spt6-set2', "replicate" = 2), "h3-set2-2" = list("barcode" = 'GTGGGATA', "group" = 'h3-set2', "ip" = 'h3', "gtype" = 'set2', "replicate" = 2), "input-wt-1" = list("barcode" = 'ATCCTATTC', "group" = 'input-wt', "ip" = 'input', "gtype" = 'wt', "replicate" = 1), "input-spt6-1" = list("barcode" = 'CGGACGTGG', "group" = 'input-spt6', "ip" = 'input', "gtype" = 'spt6', "replicate" = 1), "input-spt6-set2-1" = list("barcode" = 'GCGTTTCGA', "group" = 'input-spt6-set2', "ip" = 'input', "gtype" = 'spt6-set2', "replicate" = 1), "input-set2-1" = list("barcode" = 'TATCTCCG', "group" = 'input-set2', "ip" = 'input', "gtype" = 'set2', "replicate" = 1), "input-wt-2" = list("barcode" = 'CACAGTTGG', "group" = 'input-wt', "ip" = 'input', "gtype" = 'wt', "replicate" = 2), "input-spt6-2" = list("barcode" = 'GTGACTACA', "group" = 'input-spt6', "ip" = 'input', "gtype" = 'spt6', "replicate" = 2), "input-spt6-set2-2" = list("barcode" = 'TGAGAGTG', "group" = 'input-spt6-set2', "ip" = 'input', "gtype" = 'spt6-set2', "replicate" = 2), "input-set2-2" = list("barcode" = 'AATGCTGAC', "group" = 'input-set2', "ip" = 'input', "gtype" = 'set2', "replicate" = 2)), "fastq" = 'Undetermined_S0.R1.fastq.gz', "demultiplex" = list("mismatch" = 1), "cutadapt" = list("qual_cutoff" = 20, "min_length" = 5), "threads" = 4, "imagetype" = list("plot_spikein_percentage" = 'png', "plot_correlations" = 'png', "plot_heatmap" = 'png', "multi_metagene" = 'png'), "normalizations" = list("ha" = 'pol2', "pol2" = 'input', "h3k36me3" = 'h3', "h3k36me2" = 'h3', "h3" = 'input'), "annotations" = list("align-tss" = list("refpoint" = 'TSS', "upstream" = 300, "downstream" = 4000, "binstat" = 'mean', "binsize" = 20, "nan_afterend" = 'y', "sort" = 'ascend', "sortby" = 'region_length', "missingdatacolor" = 1, "reference_point" = 'y', "bodylength" = 0, "bed_suffix" = 'polIItranscripts-adjustedTSS'), "align-tes" = list("refpoint" = 'TES', "upstream" = 4000, "downstream" = 300, "binstat" = 'mean', "binsize" = 20, "nan_afterend" = 'y', "sort" = 'ascend', "sortby" = 'region_length', "missingdatacolor" = 1, "reference_point" = 'y', "bodylength" = 0, "bed_suffix" = 'polIItranscripts-adjustedTSS'), "align-tss-tes" = list("refpoint" = 'TSS-TES', "upstream" = 500, "downstream" = 500, "binstat" = 'mean', "binsize" = 20, "nan_afterend" = 'y', "sort" = 'ascend', "sortby" = 'region_length', "missingdatacolor" = 1, "reference_point" = 'n', "bodylength" = 1000, "bed_suffix" = 'polIItranscripts-adjustedTSS')), "callpeaks" = list("h3k36me3" = 'input', "h3k36me2" = 'input', "ha" = 'input', "pol2" = 'input', "h3" = 'input'), "macs2" = list("fdr" = 0.05, "env_name" = 'python2')),
    rule = 'multi_metagene'
)
######## Original script #########
library(ggplot2)
library(tidyverse)
library(viridis)
library(stringr)
options(scipen = 100)

ip = snakemake@wildcards[["ip"]]
filename_chip = snakemake@input[["mat"]]
window = snakemake@params[["binsize"]]
upstream = snakemake@params[["upstream"]]
downstream = snakemake@params[["downstream"]]
bodylength = snakemake@params[["bodylength"]]
align = snakemake@params[["align"]]
outpath = snakemake@output[["matname"]]
plotname = snakemake@output[["plotname"]]

#Function to get IP from the column name
get_ip = function(string)
{
  start=0;end=0;
  start = head(gregexpr("-",string)[[1]],1) + 1
  end = head(gregexpr("_",string)[[1]],1) - 1
  ip= substr(string, start, end)
  return(ip)
}
#Breadth of the matrix belonging to each sample
feature_length = (upstream + downstream + bodylength)/window

#Define where tick marks are required
if(bodylength==0)
{
  ticks = c(0,upstream/window+0.5,feature_length)
  labels = c(-upstream,align,downstream)
} else {
  ticks = c(0,upstream/window+0.5,(upstream+bodylength)/window+0.5,feature_length)
  labels = c(-upstream,"TSS","TES",downstream)
}


#Read in the matrix: The input matrix is generated by deeptools computematrix
df_chip = filename_chip %>% read_tsv(skip=3, col_names = F)
nsamples = ncol(df_chip)/feature_length

#Get column names separately
col_line = filename_chip %>% read.table(head=F, skip =2, nrows=1, stringsAsFactors=F)
col_line = col_line[-1]
names = unlist(lapply(col_line, get_ip))
names = paste(names, rep(1:feature_length,nsamples), sep="_")
colnames(df_chip) = names


#Melt the matrix
melt_matrix = gather(df_chip, key=name, value=value)

#Split the "name" column to get gtype and x value
sample_and_x = as.data.frame(melt_matrix$name %>% str_split_fixed("_",2),stringsAsFactors=F)
colnames(sample_and_x) = c("sample","x")
gtype = str_split_fixed(sample_and_x$sample,"-[:digit:]",2)[,1]
melt_matrix = cbind(melt_matrix,sample_and_x,gtype)
melt_matrix$x = as.numeric(melt_matrix$x)

print(ticks)
print(labels)

#Plot the metagene
g=ggplot(melt_matrix, aes(x,value, color = gtype, group = sample))+
  stat_summary(fun.data=mean_cl_normal, geom="smooth", alpha = 0.3)+
  ylab("") + xlab("Distance (kb)")+ ggtitle(ip)+
  scale_x_continuous(breaks=ticks, labels=labels)+
  scale_color_viridis("",option="viridis",discrete=T)+
  ylim(0,NA)+
  theme_bw()+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.text = element_text(size=15),
        plot.title = element_text(size=25,hjust=0.5),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 15))

ggsave(plot=g, filename = plotname, width = 10, height = 5, units = "in", dpi = 300)

write_tsv(melt_matrix,outpath,col_names = T)



