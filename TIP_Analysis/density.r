

# Check if ggplot2 is installed, if so, load it. If not, install and load it
if("ggplot2" %in% rownames(installed.packages())){
  library(ggplot2)
} else {
  install.packages("ggplot2")
  library(ggplot2)
}

# Import a text file with gene positions
# Column Headers: Chr, Strt (No end or Gene name required)

tips <- read.table("/Users/shelvasha/Grif16309/TEfinder_20210726224912/TE_density.bed",sep="\t",header=T)
colnames(tips)[1:2] <- c("Chr","Start")

# Orders chromosomes to appear properly in the plot
tips$Chr <- with(tips, factor(Chr, levels=paste("chr",c(1:22,"X","Y"),sep=""), ordered=TRUE))

# Cretes density plot of genes over the provided chromosomes (or scaffolds ...)
plottedTips <- ggplot(tips) + geom_histogram(aes(x=Start),binwidth=1000000) + facet_wrap(~Chr,ncol=2) + ggtitle("TIPs density over S. bicolor genome") + xlab("Genomic position (bins 1 Mb)") + ylab("Number of TIPs")
plottedTips

# Save it to an image
png("genes.png",width=1000,height=1500)
print(plottedTips)
