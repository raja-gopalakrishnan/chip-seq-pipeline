library(dplyr)
library(ggplot2)
library(tidyverse)
library(pals)

window = snakemake@params[["binsize"]]
upstream = snakemake@params[["upstream"]]
downstream = snakemake@params[["downstream"]]
bodylength = snakemake@params[["bodylength"]]
align = snakemake@params[["align"]]
annotation = snakemake@wildcards[["annotation"]]

imagetype = snakemake@params[["imagetype"]]
ip=snakemake@wildcards[["ip"]]

gtypes= snakemake@params[["gtypes"]]

feature_length = (upstream + downstream + bodylength)/window
#Define tick marks
if(bodylength==0)
{
  ticks = c(0,upstream/window+0.5,feature_length)
  labels = c(-upstream/1000,align,downstream/1000)
  xint = upstream/window+0.5
} else {
  ticks = c(0,upstream/window+0.5,(upstream+bodylength)/window+0.5,feature_length)
  labels = c(-upstream,"TSS","TES",downstream)
  xint = c(upstream/window+0.5,(upstream+bodylength)/window+0.5)
}

if(align=="TSS" && downstream >=1000)
{
  ticks = c(upstream/window,(upstream/window)+seq(1000,downstream, length.out = length(seq(1,downstream/1000)))/window) +0.5
  labels = c(align,seq(1,downstream/1000))
}

for(i in gtypes)
{
  print(i)
  #Read in IP matrix (ouptut from deeptools)
  df1=read.table(snakemake@input[["dtfile"]], skip=3, stringsAsFactors = F)
  
  df=cbind.data.frame(nrow(df1):1,df1)
  colnames(df) = c("y",c(1:feature_length))
  
  print(df[1:5,1:10])
  #Melt matrix
  print("step 1")
  melt_final = df %>% gather(key=x, value=ip_value, -y)
  print("step 2")
  melt_final$x = as.numeric(melt_final$x)
  print("step 3")
  melt_final = melt_final %>% mutate(gtype = i)
  
  if(i==gtypes[1]){
    final_mat = melt_final
  } else {
    final_mat = rbind.data.frame(final_mat,melt_final)
  }
}

final_mat$gtype = factor(final_mat$gtype, levels = unique(final_mat$gtype))
mid = median(final_mat$ip_value, na.rm = T)
print(mid)

g=ggplot(final_mat, aes(x,y))+ 
  geom_raster(aes(fill=ip_value))+
  scale_fill_gradientn(colors=coolwarm(100), na.value="white", 
                       guide=guide_colorbar(paste0(""),
                      title.position = "bottom",barwidth = unit(2,"inches")))+
  facet_grid(.~gtype)+
  theme_bw()+
  scale_x_continuous(expand=c(0,0),breaks=ticks, labels=labels)+
  scale_y_continuous(expand=c(0,0))+
  xlab("Distance (kb)") + ylab("Genes") + ggtitle(paste(ip," occupancy\n(",nrow(df1)," genes)",sep=""))+
  theme(strip.background = element_blank(),
        axis.text.y = element_blank(),
        axis.text= element_text(family="Arial", size=15, colour = "black"),
        axis.title= element_text(family="Arial", size=15, colour = "black"),
        plot.title= element_text(size=15,family="Arial", face="bold", hjust=0.5),
        legend.direction = "horizontal", legend.position = "bottom",
        legend.title = element_text(size=15,family="Arial", face="bold"),
        legend.text = element_text(family="Arial", size=15),
        strip.text.x = element_text(family="Arial", size=15),
        panel.spacing.x=unit(0.3,"inches"),
        axis.ticks.x=element_line(size=1),
        axis.ticks.length = unit(0.1,"inches"),
        axis.ticks.y=element_blank(),
        legend.title.align=0.5)

if(length(gtypes)==1)
{
  ggsave(filename = snakemake@output[["plot"]], plot = g, width = unit(3,"inches"), height = unit(5.6, "inches"),dpi=300)
} else {
  ggsave(filename = snakemake@output[["plot"]], plot = g, width = unit(2*length(gtypes),"inches"), height = unit(5.6, "inches"),dpi=300)
}
