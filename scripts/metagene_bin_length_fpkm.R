library(tidyverse)
library(viridis)
library(tibble)
library(dplyr)

file_ip = snakemake@input[["mat_ip"]]
window = snakemake@params[["binsize"]]
upstream = snakemake@params[["upstream"]]
downstream = snakemake@params[["downstream"]]
bodylength = snakemake@params[["bodylength"]]
align = snakemake@params[["align"]]
sample = snakemake@wildcards[["sample"]]
imagetype = snakemake@params[["imagetype"]]

length_outdir = snakemake@params[["length_outdir"]]
fpkm_outdir = snakemake@params[["fpkm_outdir"]]

gtype=snakemake@params[["gtype"]]
replicate=snakemake@params[["replicate"]]
ip=snakemake@params[["ip"]]
length_bins=snakemake@params[["length_bins"]]
fpkm_bins=snakemake@params[["fpkm_bins"]]

feature_length = (upstream + downstream + bodylength)/window

dir.create(length_outdir, recursive=T)
dir.create(fpkm_outdir, recursive=T)
#Read in IP matrix
df = read.table(file_ip, skip=1, stringsAsFactors = F)
#Read in FPKM table
fpkm=read.table(snakemake@params[["fpkm_table"]], head=T, stringsAsFactors = F)

print("Input matrix read. FPKM table read")

#Define tick marks
if(bodylength==0)
{
  ticks = c(0,upstream/window+0.5,feature_length)
  labels = c(-upstream,align,downstream)
  xint = upstream/window+0.5
} else {
  ticks = c(0,upstream/window+0.5,(upstream+bodylength)/window+0.5,feature_length)
  labels = c(-upstream,"TSS","TES",downstream)
  xint = c(upstream/window+0.5,(upstream+bodylength)/window+0.5)
}

#Find fpkm values for every genes from the "fpkm" file and add it as a new column to the df data frame
gene_match=NULL
for(i in df[,4])
{
  gene_match = c(gene_match,match(i,fpkm$common_name))
}
df=add_column(df,((fpkm[gene_match,6]+fpkm[gene_match,7])/2),.after = 6)

#Calculate length of genes and bind into n quantiles as defined by the groups variable
df2=cbind.data.frame((df[3]-df[2]),df[-1:-6])
colnames(df2)=c("length","fpkm",c(1:feature_length))
df2 = df2 %>% mutate(length_quantile=ntile(length,length_bins), fpkm_quantile = ntile(fpkm,fpkm_bins))

#Melt the matrix
melt_matrix = gather(df2, key=name, value=value, -length, -fpkm, -length_quantile, -fpkm_quantile)
melt_matrix$name = as.numeric(melt_matrix$name)

print("IP matrix melted")

divplot = function(control)
{
  print(control)
  #Read in control matrix
  file_control = paste(snakemake@params[["prefix"]],control,"-",gtype,"-",replicate,snakemake@params[["suffix"]],sep="")
  length_outpath = paste(snakemake@params[["length_outpath"]],control,"_length_bins.",imagetype,sep="")
  fpkm_outpath = paste(snakemake@params[["fpkm_outpath"]],control,"_fpkm_bins.",imagetype,sep="")
  
  df_control = read.table(file_control, skip=1, stringsAsFactors = F)
  
  #Do the same calculation of length of genes and melting for the control matrix
  df3=cbind.data.frame((df_control[3]-df_control[2]),df_control[-1:-6])
  colnames(df3)=c("length",c(1:feature_length))
  
  melt_control = gather(df3, key=name, value=value, -length)
  melt_control$name = as.numeric(melt_control$name)
  
  #Combine the IP and control matrices into one matrix (add an extra column for the control matrix and filter out the NAs)
  melt_final = cbind.data.frame(melt_matrix,melt_control$value)
  colnames(melt_final)[7] = "value2"
  melt_final = melt_final %>% filter(!is.na(value)) %>% filter(!is.na(value2)) %>% filter(!is.na(fpkm_quantile))
  
  #Give a pseudocount of 1 and obtain the log2(IP/control ratio): Normalization step
  new_mat = melt_final %>% mutate(div = log2(value+1)-log2(value2+1))
  
  #Legend labels for the FPKM plot
  fpkm_labels = NULL
  for(i in 1:fpkm_bins)
  {
    fpkm_min = round(min(new_mat$fpkm[which(new_mat$fpkm_quantile==i)]),1)
    fpkm_max = round(max(new_mat$fpkm[which(new_mat$fpkm_quantile==i)]),1)
    fpkm_labels = c(fpkm_labels, paste(fpkm_min," - ",fpkm_max," (",sum(df2$fpkm_quantile==i, na.rm = T)," genes)",sep=""))
  }
  
  #Legend labels for the length binned plot
  length_labels = NULL
  for(i in 1:length_bins)
  {
    length_min = round(min(new_mat$length[which(new_mat$length_quantile==i)])/1000,1)
    length_max = round(max(new_mat$length[which(new_mat$length_quantile==i)])/1000,1)
    length_labels = c(length_labels, paste(length_min," - ",length_max," kb (",sum(df2$length_quantile==i, na.rm = T)," genes)",sep=""))
  }
  
  #Make the plot for the normalized datasets: grouping by FPKM
  g1=ggplot(new_mat, aes(name,div, color= as.factor(fpkm_quantile)))+
    stat_summary(fun.data=mean_cl_normal, geom="smooth", alpha = 0.3)+
    ylab(bquote('log'[2]~.(paste(ip,"/",control,sep="")))) + xlab("Distance (kb)")+ ggtitle(sample)+
    scale_x_continuous(breaks=ticks, labels=labels)+
    scale_color_viridis("FPKM",option="viridis",discrete=T,labels=fpkm_labels)+
    theme_bw()+
    theme(panel.grid.minor=element_blank(),
          axis.text = element_text(size=15),
          plot.title = element_text(size=25,hjust=0.5),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15))
  
  #Make the plot for the normalized datasets: grouping by FPKM
  g2=ggplot(new_mat, aes(name,div, color= as.factor(length_quantile)))+
    stat_summary(fun.data=mean_cl_normal, geom="smooth", alpha = 0.3)+
    ylab(bquote('log'[2]~.(paste(ip,"/",control,sep="")))) + xlab("Distance (kb)")+ ggtitle(sample)+
    scale_x_continuous(breaks=ticks, labels=labels)+
    scale_color_viridis("Gene length",option="viridis",discrete=T,labels=length_labels)+
    theme_bw()+
    theme(panel.grid.minor=element_blank(),
          axis.text = element_text(size=15),
          plot.title = element_text(size=25,hjust=0.5),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15))
  
  #Save the plots
  ggsave(filename = fpkm_outpath, plot = g1, width = 10, height = 5, units = "in")
  ggsave(filename = length_outpath, plot = g2, width = 10, height = 5, units = "in")
}

for(i in snakemake@params[["control"]]) { divplot(i) }

write.table("Done!",snakemake@output[[1]])
