library(ggplot2)
library(reshape2)
library(dplyr)

tab = NULL
columns = c("S.cerevisiae","S.pombe","total","sample","ip")
#Import all text files containing read numbers and paste them together to make a file
for (i in 1:length(snakemake@input[["files"]]))
{
	temp = t(read.table(snakemake@input[["files"]][i], stringsAsFactors=F, header =F))
	temp = cbind.data.frame(temp,names(snakemake@params[["samples"]])[i],snakemake@params[["ip"]][i])
	colnames(temp)=columns
	tab = rbind.data.frame(tab, temp)
}
#Write the table with read values for all samples 
write.table(tab,snakemake@output[["all_txt"]], quote=F, row.names=F, sep="\t")

#Calculate proportion of cerevisiae and spike in reads for each library by dividing by total number of reads
tab[,1:3] = tab[,1:3]/tab[,3]

#Create a melted table amenable for ggplot (Exclude column 3 which is total number of reads)
melt_tab = melt(tab[,-3])

unique_ips=unique(snakemake@params[["ip"]])
#Write a function to make the plot. Index will vary from 1 to 6 and make 6 plots for all different IPs and input samples
plot_prop = function(index)
{
	#Filter the data to just have the IP samples required
	IP= unique_ips[index]
	z = filter(melt_tab, ip==IP)
	#Make the plot
	g = ggplot(z,aes(x=sample, y=value, fill=variable)) +geom_bar(stat = "identity")+theme_bw()+
		  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
		        axis.title.x = element_text(size = 50), axis.title.y = element_blank(),
		        axis.text.x = element_text(size = 60), axis.text.y = element_text(size=60),
		        legend.text = element_text(size=50),legend.key.size = unit(7,"lines"),
		        legend.title = element_blank())+
		  scale_fill_manual(values=c(rgb(252/255,176/255,64/255),rgb(39/255,34/255,98/255)))+
		  scale_y_continuous(expand=c(0,0), breaks = c(0.25,0.5,0.75,1))+
		  ylab("Proportion")+coord_flip()

	#Save the plot
	ggsave(filename = snakemake@output[["ip_plots"]][index], plot=g,width=8.3, height=5, units="in")
}

#Loop through to make plots for all ips
for(i in 1:length(unique_ips)) { plot_prop(i) }



