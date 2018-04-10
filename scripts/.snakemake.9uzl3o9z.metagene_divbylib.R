
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
    input = list('matrix/group_ip_melted/align-tss/Scer_normRPM/ha_Scer_normRPM_align-tss_meltedmat.txt', 'matrix/group_ip_melted/align-tss/Spom_normRPM/ha_Spom_normRPM_align-tss_meltedmat.txt', 'matrix/group_ip_melted/align-tes/Scer_normRPM/ha_Scer_normRPM_align-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tes/Spom_normRPM/ha_Spom_normRPM_align-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss-tes/Scer_normRPM/ha_Scer_normRPM_align-tss-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss-tes/Spom_normRPM/ha_Spom_normRPM_align-tss-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss/Scer_normRPM/pol2_Scer_normRPM_align-tss_meltedmat.txt', 'matrix/group_ip_melted/align-tss/Spom_normRPM/pol2_Spom_normRPM_align-tss_meltedmat.txt', 'matrix/group_ip_melted/align-tes/Scer_normRPM/pol2_Scer_normRPM_align-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tes/Spom_normRPM/pol2_Spom_normRPM_align-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss-tes/Scer_normRPM/pol2_Scer_normRPM_align-tss-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss-tes/Spom_normRPM/pol2_Spom_normRPM_align-tss-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss/Scer_normRPM/h3k36me2_Scer_normRPM_align-tss_meltedmat.txt', 'matrix/group_ip_melted/align-tss/Spom_normRPM/h3k36me2_Spom_normRPM_align-tss_meltedmat.txt', 'matrix/group_ip_melted/align-tes/Scer_normRPM/h3k36me2_Scer_normRPM_align-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tes/Spom_normRPM/h3k36me2_Spom_normRPM_align-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss-tes/Scer_normRPM/h3k36me2_Scer_normRPM_align-tss-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss-tes/Spom_normRPM/h3k36me2_Spom_normRPM_align-tss-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss/Scer_normRPM/h3k36me3_Scer_normRPM_align-tss_meltedmat.txt', 'matrix/group_ip_melted/align-tss/Spom_normRPM/h3k36me3_Spom_normRPM_align-tss_meltedmat.txt', 'matrix/group_ip_melted/align-tes/Scer_normRPM/h3k36me3_Scer_normRPM_align-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tes/Spom_normRPM/h3k36me3_Spom_normRPM_align-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss-tes/Scer_normRPM/h3k36me3_Scer_normRPM_align-tss-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss-tes/Spom_normRPM/h3k36me3_Spom_normRPM_align-tss-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss/Scer_normRPM/h3_Scer_normRPM_align-tss_meltedmat.txt', 'matrix/group_ip_melted/align-tss/Spom_normRPM/h3_Spom_normRPM_align-tss_meltedmat.txt', 'matrix/group_ip_melted/align-tes/Scer_normRPM/h3_Scer_normRPM_align-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tes/Spom_normRPM/h3_Spom_normRPM_align-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss-tes/Scer_normRPM/h3_Scer_normRPM_align-tss-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss-tes/Spom_normRPM/h3_Spom_normRPM_align-tss-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss/Scer_normRPM/input_Scer_normRPM_align-tss_meltedmat.txt', 'matrix/group_ip_melted/align-tss/Spom_normRPM/input_Spom_normRPM_align-tss_meltedmat.txt', 'matrix/group_ip_melted/align-tes/Scer_normRPM/input_Scer_normRPM_align-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tes/Spom_normRPM/input_Spom_normRPM_align-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss-tes/Scer_normRPM/input_Scer_normRPM_align-tss-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss-tes/Spom_normRPM/input_Spom_normRPM_align-tss-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss/Scer_normSI/ha_Scer_normSI_align-tss_meltedmat.txt', 'matrix/group_ip_melted/align-tss/Spom_normSI/ha_Spom_normSI_align-tss_meltedmat.txt', 'matrix/group_ip_melted/align-tes/Scer_normSI/ha_Scer_normSI_align-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tes/Spom_normSI/ha_Spom_normSI_align-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss-tes/Scer_normSI/ha_Scer_normSI_align-tss-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss-tes/Spom_normSI/ha_Spom_normSI_align-tss-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss/Scer_normSI/pol2_Scer_normSI_align-tss_meltedmat.txt', 'matrix/group_ip_melted/align-tss/Spom_normSI/pol2_Spom_normSI_align-tss_meltedmat.txt', 'matrix/group_ip_melted/align-tes/Scer_normSI/pol2_Scer_normSI_align-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tes/Spom_normSI/pol2_Spom_normSI_align-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss-tes/Scer_normSI/pol2_Scer_normSI_align-tss-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss-tes/Spom_normSI/pol2_Spom_normSI_align-tss-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss/Scer_normSI/h3k36me2_Scer_normSI_align-tss_meltedmat.txt', 'matrix/group_ip_melted/align-tss/Spom_normSI/h3k36me2_Spom_normSI_align-tss_meltedmat.txt', 'matrix/group_ip_melted/align-tes/Scer_normSI/h3k36me2_Scer_normSI_align-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tes/Spom_normSI/h3k36me2_Spom_normSI_align-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss-tes/Scer_normSI/h3k36me2_Scer_normSI_align-tss-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss-tes/Spom_normSI/h3k36me2_Spom_normSI_align-tss-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss/Scer_normSI/h3k36me3_Scer_normSI_align-tss_meltedmat.txt', 'matrix/group_ip_melted/align-tss/Spom_normSI/h3k36me3_Spom_normSI_align-tss_meltedmat.txt', 'matrix/group_ip_melted/align-tes/Scer_normSI/h3k36me3_Scer_normSI_align-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tes/Spom_normSI/h3k36me3_Spom_normSI_align-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss-tes/Scer_normSI/h3k36me3_Scer_normSI_align-tss-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss-tes/Spom_normSI/h3k36me3_Spom_normSI_align-tss-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss/Scer_normSI/h3_Scer_normSI_align-tss_meltedmat.txt', 'matrix/group_ip_melted/align-tss/Spom_normSI/h3_Spom_normSI_align-tss_meltedmat.txt', 'matrix/group_ip_melted/align-tes/Scer_normSI/h3_Scer_normSI_align-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tes/Spom_normSI/h3_Spom_normSI_align-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss-tes/Scer_normSI/h3_Scer_normSI_align-tss-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss-tes/Spom_normSI/h3_Spom_normSI_align-tss-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss/Scer_normSI/input_Scer_normSI_align-tss_meltedmat.txt', 'matrix/group_ip_melted/align-tss/Spom_normSI/input_Spom_normSI_align-tss_meltedmat.txt', 'matrix/group_ip_melted/align-tes/Scer_normSI/input_Scer_normSI_align-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tes/Spom_normSI/input_Spom_normSI_align-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss-tes/Scer_normSI/input_Scer_normSI_align-tss-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss-tes/Spom_normSI/input_Spom_normSI_align-tss-tes_meltedmat.txt', 'matrix/group_ip_melted/align-tss/Spom_normRPM/h3_Spom_normRPM_align-tss_meltedmat.txt', "mat_ip" = 'matrix/group_ip_melted/align-tss/Spom_normRPM/h3_Spom_normRPM_align-tss_meltedmat.txt'),
    output = list('plots/metagene/group_ip_div_lib/align-tss/Spom_normRPM/h3_Spom_normRPM_align-tss_div___snakemake_dynamic___metagene.png'),
    params = list('matrix/group_ip_melted/align-tss/Spom_normRPM/', '_Spom_normRPM_align-tss_meltedmat.txt', c('input'), 20, 300, 4000, 0, 'TSS', 'spt6', 'plots/metagene/group_ip_div_lib/align-tss/Spom_normRPM/h3_Spom_normRPM_align-tss_div_', 'png', "prefix" = 'matrix/group_ip_melted/align-tss/Spom_normRPM/', "suffix" = '_Spom_normRPM_align-tss_meltedmat.txt', "control" = c('input'), "binsize" = 20, "upstream" = 300, "downstream" = 4000, "bodylength" = 0, "align" = 'TSS', "exclude" = 'spt6', "outpath" = 'plots/metagene/group_ip_div_lib/align-tss/Spom_normRPM/h3_Spom_normRPM_align-tss_div_', "imagetype" = 'png'),
    wildcards = list('align-tss', 'Spom', 'RPM', 'h3', '__snakemake_dynamic__', "annotation" = 'align-tss', "species" = 'Spom', "norm" = 'RPM', "ip" = 'h3', "control" = '__snakemake_dynamic__'),
    threads = 1,
    log = list(),
    resources = list(),
    config = list("samples" = list("ha-wt-1" = list("barcode" = 'ATCACG', "group" = 'ha-wt', "ip" = 'ha', "gtype" = 'wt', "replicate" = 1), "ha-spt6-1" = list("barcode" = 'CGATGT', "group" = 'ha-spt6', "ip" = 'ha', "gtype" = 'spt6', "replicate" = 1), "ha-spt6-set2-1" = list("barcode" = 'TTAGGC', "group" = 'ha-spt6-set2', "ip" = 'ha', "gtype" = 'spt6-set2', "replicate" = 1), "ha-set2-1" = list("barcode" = 'TGACCA', "group" = 'ha-set2', "ip" = 'ha', "gtype" = 'set2', "replicate" = 1), "ha-wt-2" = list("barcode" = 'ACAGTG', "group" = 'ha-wt', "ip" = 'ha', "gtype" = 'wt', "replicate" = 2), "ha-spt6-2" = list("barcode" = 'GCCAAT', "group" = 'ha-spt6', "ip" = 'ha', "gtype" = 'spt6', "replicate" = 2), "ha-spt6-set2-2" = list("barcode" = 'CAGATC', "group" = 'ha-spt6-set2', "ip" = 'ha', "gtype" = 'spt6-set2', "replicate" = 2), "ha-set2-2" = list("barcode" = 'ACTTGA', "group" = 'ha-set2', "ip" = 'ha', "gtype" = 'set2', "replicate" = 2), "pol2-wt-1" = list("barcode" = 'GATCAG', "group" = 'pol2-wt', "ip" = 'pol2', "gtype" = 'wt', "replicate" = 1), "pol2-spt6-1" = list("barcode" = 'TAGCTT', "group" = 'pol2-spt6', "ip" = 'pol2', "gtype" = 'spt6', "replicate" = 1), "pol2-spt6-set2-1" = list("barcode" = 'GGCTAC', "group" = 'pol2-spt6-set2', "ip" = 'pol2', "gtype" = 'spt6-set2', "replicate" = 1), "pol2-set2-1" = list("barcode" = 'CTTGTA', "group" = 'pol2-set2', "ip" = 'pol2', "gtype" = 'set2', "replicate" = 1), "pol2-wt-2" = list("barcode" = 'ATATAGGA', "group" = 'pol2-wt', "ip" = 'pol2', "gtype" = 'wt', "replicate" = 2), "pol2-spt6-2" = list("barcode" = 'AACCGTGT', "group" = 'pol2-spt6', "ip" = 'pol2', "gtype" = 'spt6', "replicate" = 2), "pol2-spt6-set2-2" = list("barcode" = 'AGGTCAGT', "group" = 'pol2-spt6-set2', "ip" = 'pol2', "gtype" = 'spt6-set2', "replicate" = 2), "pol2-set2-2" = list("barcode" = 'CTCTGTCT', "group" = 'pol2-set2', "ip" = 'pol2', "gtype" = 'set2', "replicate" = 2), "h3k36me2-wt-1" = list("barcode" = 'CCATACAC', "group" = 'h3k36me2-wt', "ip" = 'h3k36me2', "gtype" = 'wt', "replicate" = 1), "h3k36me2-spt6-1" = list("barcode" = 'CGCATTAA', "group" = 'h3k36me2-spt6', "ip" = 'h3k36me2', "gtype" = 'spt6', "replicate" = 1), "h3k36me2-spt6-set2-1" = list("barcode" = 'GTCTACAT', "group" = 'h3k36me2-spt6-set2', "ip" = 'h3k36me2', "gtype" = 'spt6-set2', "replicate" = 1), "h3k36me2-set2-1" = list("barcode" = 'GAGTTAAC', "group" = 'h3k36me2-set2', "ip" = 'h3k36me2', "gtype" = 'set2', "replicate" = 1), "h3k36me2-wt-2" = list("barcode" = 'GCAGCCTC', "group" = 'h3k36me2-wt', "ip" = 'h3k36me2', "gtype" = 'wt', "replicate" = 2), "h3k36me2-spt6-2" = list("barcode" = 'TCGCGTAC', "group" = 'h3k36me2-spt6', "ip" = 'h3k36me2', "gtype" = 'spt6', "replicate" = 2), "h3k36me2-spt6-set2-2" = list("barcode" = 'TATACCGT', "group" = 'h3k36me2-spt6-set2', "ip" = 'h3k36me2', "gtype" = 'spt6-set2', "replicate" = 2), "h3k36me2-set2-2" = list("barcode" = 'TGCGGTTA', "group" = 'h3k36me2-set2', "ip" = 'h3k36me2', "gtype" = 'set2', "replicate" = 2), "h3k36me3-wt-1" = list("barcode" = 'AACACCTAC', "group" = 'h3k36me3-wt', "ip" = 'h3k36me3', "gtype" = 'wt', "replicate" = 1), "h3k36me3-spt6-1" = list("barcode" = 'CCTTTACAG', "group" = 'h3k36me3-spt6', "ip" = 'h3k36me3', "gtype" = 'spt6', "replicate" = 1), "h3k36me3-spt6-set2-1" = list("barcode" = 'GGTCCTTGA', "group" = 'h3k36me3-spt6-set2', "ip" = 'h3k36me3', "gtype" = 'spt6-set2', "replicate" = 1), "h3k36me3-set2-1" = list("barcode" = 'TTGAGTGT', "group" = 'h3k36me3-set2', "ip" = 'h3k36me3', "gtype" = 'set2', "replicate" = 1), "h3k36me3-wt-2" = list("barcode" = 'ACTAACTGC', "group" = 'h3k36me3-wt', "ip" = 'h3k36me3', "gtype" = 'wt', "replicate" = 2), "h3k36me3-spt6-2" = list("barcode" = 'CAGGAGGCG', "group" = 'h3k36me3-spt6', "ip" = 'h3k36me3', "gtype" = 'spt6', "replicate" = 2), "h3k36me3-spt6-set2-2" = list("barcode" = 'GTTGTCCCA', "group" = 'h3k36me3-spt6-set2', "ip" = 'h3k36me3', "gtype" = 'spt6-set2', "replicate" = 2), "h3k36me3-set2-2" = list("barcode" = 'TGACGCAT', "group" = 'h3k36me3-set2', "ip" = 'h3k36me3', "gtype" = 'set2', "replicate" = 2), "h3-wt-1" = list("barcode" = 'ATCGCCAGC', "group" = 'h3-wt', "ip" = 'h3', "gtype" = 'wt', "replicate" = 1), "h3-spt6-1" = list("barcode" = 'CATTCCAAG', "group" = 'h3-spt6', "ip" = 'h3', "gtype" = 'spt6', "replicate" = 1), "h3-spt6-set2-1" = list("barcode" = 'GCAAGTAGA', "group" = 'h3-spt6-set2', "ip" = 'h3', "gtype" = 'spt6-set2', "replicate" = 1), "h3-set2-1" = list("barcode" = 'TGATCCGA', "group" = 'h3-set2', "ip" = 'h3', "gtype" = 'set2', "replicate" = 1), "h3-wt-2" = list("barcode" = 'ACGTAGCTC', "group" = 'h3-wt', "ip" = 'h3', "gtype" = 'wt', "replicate" = 2), "h3-spt6-2" = list("barcode" = 'CGAACTGTG', "group" = 'h3-spt6', "ip" = 'h3', "gtype" = 'spt6', "replicate" = 2), "h3-spt6-set2-2" = list("barcode" = 'TAGCTAGTA', "group" = 'h3-spt6-set2', "ip" = 'h3', "gtype" = 'spt6-set2', "replicate" = 2), "h3-set2-2" = list("barcode" = 'GTGGGATA', "group" = 'h3-set2', "ip" = 'h3', "gtype" = 'set2', "replicate" = 2), "input-wt-1" = list("barcode" = 'ATCCTATTC', "group" = 'input-wt', "ip" = 'input', "gtype" = 'wt', "replicate" = 1), "input-spt6-1" = list("barcode" = 'CGGACGTGG', "group" = 'input-spt6', "ip" = 'input', "gtype" = 'spt6', "replicate" = 1), "input-spt6-set2-1" = list("barcode" = 'GCGTTTCGA', "group" = 'input-spt6-set2', "ip" = 'input', "gtype" = 'spt6-set2', "replicate" = 1), "input-set2-1" = list("barcode" = 'TATCTCCG', "group" = 'input-set2', "ip" = 'input', "gtype" = 'set2', "replicate" = 1), "input-wt-2" = list("barcode" = 'CACAGTTGG', "group" = 'input-wt', "ip" = 'input', "gtype" = 'wt', "replicate" = 2), "input-spt6-2" = list("barcode" = 'GTGACTACA', "group" = 'input-spt6', "ip" = 'input', "gtype" = 'spt6', "replicate" = 2), "input-spt6-set2-2" = list("barcode" = 'TGAGAGTG', "group" = 'input-spt6-set2', "ip" = 'input', "gtype" = 'spt6-set2', "replicate" = 2), "input-set2-2" = list("barcode" = 'AATGCTGAC', "group" = 'input-set2', "ip" = 'input', "gtype" = 'set2', "replicate" = 2)), "fastq" = 'Undetermined_S0.R1.fastq.gz', "demultiplex" = list("mismatch" = 1), "cutadapt" = list("qual_cutoff" = 20, "min_length" = 5), "threads" = 4, "imagetype" = list("plot_spikein_percentage" = 'png', "plot_correlations" = 'png', "plot_heatmap" = 'png', "multi_metagene" = 'png', "multi_metagene_divbylib" = 'png'), "exclude" = list("multi_metagene_divbylib" = 'spt6'), "normalizations" = list("ha" = c('input', 'pol2'), "pol2" = c('input'), "h3k36me3" = c('input', 'h3', 'set2'), "h3k36me2" = c('input', 'h3', 'set2'), "h3" = c('input')), "annotations" = list("align-tss" = list("refpoint" = 'TSS', "upstream" = 300, "downstream" = 4000, "binstat" = 'mean', "binsize" = 20, "nan_afterend" = 'y', "sort" = 'ascend', "sortby" = 'region_length', "missingdatacolor" = 1, "reference_point" = 'y', "bodylength" = 0, "bed_suffix" = 'polIItranscripts-adjustedTSS'), "align-tes" = list("refpoint" = 'TES', "upstream" = 4000, "downstream" = 300, "binstat" = 'mean', "binsize" = 20, "nan_afterend" = 'y', "sort" = 'ascend', "sortby" = 'region_length', "missingdatacolor" = 1, "reference_point" = 'y', "bodylength" = 0, "bed_suffix" = 'polIItranscripts-adjustedTSS'), "align-tss-tes" = list("refpoint" = 'TSS-TES', "upstream" = 500, "downstream" = 500, "binstat" = 'mean', "binsize" = 20, "nan_afterend" = 'y', "sort" = 'ascend', "sortby" = 'region_length', "missingdatacolor" = 1, "reference_point" = 'n', "bodylength" = 1000, "bed_suffix" = 'polIItranscripts-adjustedTSS')), "callpeaks" = list("h3k36me3" = 'input', "h3k36me2" = 'input', "ha" = 'input', "pol2" = 'input', "h3" = 'input'), "macs2" = list("fdr" = 0.05, "env_name" = 'python2')),
    rule = 'mutli_metagene_divbylib'
)
######## Original script #########
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

divplot = function(control)
{
  file_control = paste(snakemake@params[["prefix"]],control,snakemake@params[["suffix"]],sep="")
  outpath = paste(snakemake@params[["outpath"]],control,"_metagene.",imagetype,sep="")
  outpath_exclude = paste(snakemake@params[["outpath"]],control,"_no_",exclude,"_metagene.",imagetype,sep="")
  
  melt_matrix_control = read_tsv(file_control)
  
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
  }  
}

for(i in snakemake@params[["control"]]) { divplot(i) }