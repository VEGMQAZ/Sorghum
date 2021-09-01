# Check if ggplot2 is installed, if so, load it. If not, install and load it

if("ggplot2" %in% rownames(installed.packages())){
  library(ggplot2)
} else {
  install.packages("ggplot2")
  library(ggplot2)
}

# Import a text file with gene positions
# Column Headers: Chr, Strt (No end or Gene name required)

genes <- read.table("Grif16309_density.bed",sep="\t",header=T)

# make sure the chromosomes are ordered in the way you want them to appear in the plot

genes$chr <- with(genes, factor(chr, levels=paste("chr",c(1:22,"X","Y"),sep=""), ordered=TRUE))

# make a density plot of genes over the provided chromosomes (or scaffolds ...)

plottedGenes <- ggplot(genes) + geom_histogram(aes(x=pos),binwidth=1000000) + facet_wrap(~chr,ncol=2) + ggtitle("RefSeq genes density over human genome 19") + xlab("Genomic position (bins 1 Mb)") + ylab("Number of genes")

# Save it to an image
png("genes.png",width=1000,height=1500)
print(plottedGenes)
dev.off()