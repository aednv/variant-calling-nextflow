#!/usr/bin/env Rscript

library("tidyverse")
library("scales")

args <- commandArgs(TRUE)
combinedVcf <- as.character(args[1])

chr_name <- scan("chr_names.txt", what="character")
chr_length <- scan("chr_lengths.txt", what="integer")
df <- read_tsv("variantLocations.tsv", trim_ws=TRUE)
df$CHROM <- as.factor(df$CHROM)
df$POS <- as.integer(df$POS)
plot <- ggplot(df, aes(x=POS)) + geom_histogram(binwidth=500000) + scale_x_continuous(labels = comma, breaks = round(seq(min(0), max(df$POS), by=10000000), 5000000)) + 
	facet_wrap(~CHROM) + theme_bw() + theme( axis.text = element_text( angle = 45, size = 4)) + labs(title = paste(combinedVcf, " homozygous variants"), y = "count (binwidth = 500,000 bp)", x = "chromosome position (bp)")
ggsave(paste0(combinedVcf, ".variant.plot.all.pdf"), device="pdf")
for (i in 1:length(chr_name)) {
	df_chr <- df %>% filter(CHROM == chr_name[i])
	if (is.finite(max(df_chr$POS)) == FALSE) {
		break
	}
	chr_plot <- ggplot(df_chr, aes(x=POS)) + geom_histogram(binwidth=500000) + scale_x_continuous(labels = comma, breaks = round(seq(min(0), max(df_chr$POS), by=10000000), as.numeric(chr_length[i]))) + theme_bw() + 
		labs(title = paste(combinedVcf, chr_name[i], sep=" "), y = "count (binwidth = 500,000 bp)", x = "chromosome position (bp)")
	save_name_p1 <- paste( chr_name[i], combinedVcf, sep=".")
	save_name <- paste0( save_name_p1, ".variant.frequency.plot.pdf")
	ggsave(save_name, plot = chr_plot, device = "pdf")
}

