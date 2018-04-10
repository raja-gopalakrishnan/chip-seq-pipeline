library(tidyverse)
library(ggplot2)
library(viridis)
options(scipen=100)

file_ip = snakemake@input[["mat_ip"]]
window = snakemake@params[["binsize"]]
upstream = snakemake@params[["upstream"]]
downstream = snakemake@params[["downstream"]]
bodylength = snakemake@params[["bodylength"]]
align = snakemake@params[["align"]]
ip = snakemake@wildcards[["ip"]]
exclude = snakemake@params[["exclude"]]
imagetype = snakemake@params[["imagetype"]]
outdir = snakemake@params[["outdir"]]

dir.create(outdir, recursive=T)
feature_length = (upstream + downstream + bodylength)/window
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
melt_matrix_ip = read_tsv(file_ip)
print("IP mat read")
divplot = function(control)
{
  print(control)
  file_control = paste(snakemake@params[["prefix"]],control,snakemake@params[["suffix"]],sep="")
  outpath = paste(snakemake@params[["outpath"]],control,"_metagene.",imagetype,sep="")
  outpath_exclude = paste(snakemake@params[["outpath"]],control,"_metagene_no_",exclude,".",imagetype,sep="")
  
  melt_matrix_control = read_tsv(file_control)
  print("Control mat read")
  
  #Paste the two matrices side by side
  new_mat = cbind(melt_matrix_ip,melt_matrix_control)
  #Change column names to prevent redundancy
  colnames(new_mat)[1:5] = paste(colnames(new_mat)[1:5],"_ip",sep="")
  colnames(new_mat)[6:10] = paste(colnames(new_mat)[6:10],"_control",sep="")
  #Filter out NAs
  new_mat = new_mat %>% filter(!is.na(value_ip)) %>% filter(!is.na(value_control))
  #Take log2 of ip/control using a pseudocount of 1
  #The if statement is just a safety measure to ensure the right rows are pasted next to one another
  if(sum(new_mat$name_ip != new_mat$name_control)==0)
  {
    new_mat = new_mat %>% mutate(div = log2(value_ip+1)-log2(value_control+1))
  }
  new_mat$gtype_ip = factor(new_mat$gtype_ip, levels = unique(new_mat$gtype_ip))
  
  #Make the plot
  g1 = ggplot(new_mat, aes(x_ip,div, color = gtype_ip, group = sample_ip))+
    stat_summary(fun.data=mean_cl_normal, geom="smooth", alpha = 0.3)+
    ylab("log2 fold change") + xlab("Distance (kb)")+ ggtitle(paste(ip,"/",control,sep=""))+
    scale_x_continuous(breaks=ticks, labels=labels)+
    scale_color_viridis("",option="viridis",discrete=T)+
    geom_vline(xintercept = xint)+
    theme_bw()+
    theme(panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(),
          axis.text = element_text(size=15),
          plot.title = element_text(size=25,hjust=0.5),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 15))

  ggsave(filename = outpath, plot =g1, width = 10, height = 5, units = "in")
  print("plot1 made")
  if(exclude!="")
  {
    g2 = ggplot(new_mat %>% filter(gtype_ip!=exclude), aes(x_ip,div, color = gtype_ip, group = sample_ip))+
      stat_summary(fun.data=mean_cl_normal, geom="smooth", alpha = 0.3)+
      ylab("log2 fold change") + xlab("Distance (kb)")+ ggtitle(paste(ip,"/",control,sep=""))+
      scale_x_continuous(breaks=ticks, labels=labels)+
      scale_color_viridis("",option="viridis",discrete=T)+
      geom_vline(xintercept = xint)+
      theme_bw()+
      theme(panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(),
            axis.text = element_text(size=15),
            plot.title = element_text(size=25,hjust=0.5),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 15))
    ggsave(filename = outpath_exclude, plot =g2, width = 10, height = 5, units = "in")
    print("plot2 made")  
  } 
  print("Plots made") 
}

for(i in snakemake@params[["control"]]) { divplot(i) }

write.table("Done!",snakemake@output[[1]])