#Script written and provided by James Chuang

library(tidyverse)
library(GGally)
library(viridis)
library(forcats)

main = function(intable, pcount, samplelist, outpath){
    df = intable %>% read_tsv() %>% gather(key=sample, value=signal, -chrom, -start, -end) %>%
            filter(sample %in% samplelist)
    df$sample = factor(df$sample, ordered=TRUE)
    df = df %>% spread(sample, signal) %>% select(-chrom, -start, -end)
    
    maxsignal = max(df) + pcount
    mincor = min(cor(df, use="complete.obs"))
    plots = list()
    
    #for each row
    for (i in 1:ncol(df)){
        #for each column
        for (j in 1:ncol(df)){
            idx = ncol(df)*(i-1)+j
            #upper right (correlation)
            if (i < j){
                c = cor((df[,i]+pcount), (df[,j]+pcount), use = "complete.obs") %>% as.numeric()
                plot = ggplot(data = tibble(x=c(0,1), y=c(0,1), corr=c)) +
                        geom_rect(aes(fill=corr), xmin=0, ymin=0, xmax=1, ymax=1) +
                        annotate("text", x=0.5, y=0.5, label=sprintf("%.2f",round(c,2)), size=10) +
                        scale_x_continuous(breaks=NULL) +
                        scale_y_continuous(breaks=NULL) +
                        scale_fill_distiller(palette="Blues", limits = c(mincor,1), direction=1)
                plots[[idx]] = plot
            }
            #top left to bot right diag (density)
            else if (i == j){
                subdf = df %>% select(i) %>% gather(sample, value)
                plot = ggplot(data = subdf, aes(x=(value+pcount))) +
                        geom_density(aes(y=..scaled..), fill="black", size=1) +
                        scale_y_continuous(breaks=c(0,.5,1)) +
                        scale_x_log10(limit = c(pcount, maxsignal)) +
                        annotate("text", x=.90*maxsignal, y=0.5, hjust=1, 
                                 label=unique(subdf$sample), size=4, fontface="bold") 
                plots[[idx]] = plot
            }
            #bottom left (scatter)
            else {
                #the filtering is a quick hack to avoid the (0,0) bin taking up
                #all of the colorspace
                subdf = df %>% select(i,j) %>% gather(xsample, xvalue, -1) %>%
                            gather(ysample, yvalue, -c(2:3)) %>%
                            filter(!(xvalue < 6*pcount & yvalue < 6*pcount))
                plot = ggplot(data = subdf, aes(x=xvalue+pcount, y=yvalue+pcount)) +
                            geom_abline(intercept = 0, slope=1, color="grey80", size=.5) +
                            stat_bin_hex(geom="point", aes(color=log10(..count..)), binwidth=c(.04,.04), size=.15, shape=16) +
                            scale_fill_viridis(option="inferno") +
                            scale_color_viridis(option="inferno") +
                            scale_x_log10(limit = c(pcount, maxsignal)) +
                            scale_y_log10(limit = c(pcount, maxsignal))
                plots[[idx]] = plot
            }
        }
    }
    
    mat = ggmatrix(plots, nrow=ncol(df), ncol=ncol(df),
                   xAxisLabels = names(df), yAxisLabels = names(df)) +
                    theme_light() +
                    theme(axis.text = element_text(size=9),
                          strip.background = element_blank(),
                          strip.text = element_text(size=12, color="black", face="bold"),
                          strip.text.y = element_text(angle=180, hjust=1),
                          strip.placement="outside",
                          strip.switch.pad.grid = unit(0, "points"),
                          strip.switch.pad.wrap = unit(0, "points"))
    w = 3+ncol(df)*4
    h = 9/16*w
    
    ggsave(filename = outpath, plot=mat, width=w, height=h, units="cm")

    print(warnings())
}    
#Make the plot
main(intable = snakemake@input[[1]],
     pcount = snakemake@params[["pcount"]],
     samplelist = snakemake@params[["samplelist"]],
     outpath = snakemake@output[[1]])


