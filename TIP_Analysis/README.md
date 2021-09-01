# TIPs Detection Workflow in Sorghum

## Pre-requisites

- Bedtools
- Bedops
- Samtools
- Picard

**Repeatmasker .out conversion to .bed and formating**

- To perform TE insertion detection we need the following files:
  1. Reference fasta
  2. Reference geneome .out file from RepeatMasker (TEs)
  3. Installed pre-req’d packages.
- Running the following commands will generate the required files for analysis:

 

```    bash
#! /usr/bin/env bash
# set -x # Verbose testing outputs

# Loads necessary modules on Slurm
module load samtools
module load bedtools2
module load picard

# Load executable dependencies
PATH=$PATH$( find /projects/cooper_research1/TIP_Analysis/bin -type d -printf ":%p" )

# Changes to working directory
cd /projects/cooper_research1/TIP_Analysis/bin/TEfinder

# Sets variable of genome
SORGUM_REF='/projects/cooper_research1/TIP_Analysis/SorghumRefAnnot/Sbicolor_454_v3.0.1.fa'
REP_OUT='/projects/cooper_research1/TIP_Analysis/RepeatAnalysis/results_Sbicolor_454_v3/repeatmasker_out/Sbicolor_454_v3.0.1.fa.out'

# Removes simple repeats from RepeatMasker .out
grep -v -iE '(Motif\:[ATGC]+\-rich)|(Motif\:\([ATGC]+\)n)' $REP_OUT > TEs.gtf
TE_REF='/projects/cooper_research1/TIP_Analysis/bin/TEfinder/TEs.gtf'

# Converts .out from RepeatMasker into .bed file.
rmsk2bed < TEs.gtf bedops --merge - > Sbicolor_454_v3.0.1.bed

# # # #
## Commented out 08/12/21
## Limits to LTR/Copia elements, truncates to the first element entry in .bed
cat Sbicolor_454_v3.0.1.bed | cut -d, -f11 | grep "LTR/Copia" | head -10 > truncated.bed
# # # #

# cat Sbicolor_454_v3.0.1.bed > truncated.bed

# Removes asterisk characters
tr -d '*' < truncated.bed > Sbicolor_truncated_rmchar.bed

# Converts coordinates from .bed to fasta
bedtools getfasta -fi $SORGUM_REF -name -bed Sbicolor_truncated_rmchar.bed > Sbicolor_454_v3.0.1_TE_truncated.fa

# Removes special characters (::) from bedtools-generated fasta
sed 's/::.*//' Sbicolor_454_v3.0.1_TE_truncated.fa > Sbicolor_454_v3.0.1_TE_truncated_format.fa

# Outputs fasta entries into .txt file
grep -e ">" Sbicolor_454_v3.0.1_TE_truncated_format.fa | awk 'sub(/^>/, "")' >> TE_list.txt
TE_LIST='/projects/cooper_research1/TIP_Analysis/bin/TEfinder/TE_list.txt'


## Cleanup of intermediate files
rm -r truncated.bed
rm -r Sbicolor_truncated_rmchar.bed
rm -r Sbicolor_454_v3.0.1_TE_truncated_format.fa
rm -r Sbicolor_454_v3.0.1_TE_truncated.fa
        
```

### Detecting TE insertions

Once the setup and gathering of files is complete, we can proceed with the discovery of non-reference TE insertions using the following script:

```bash
## Running TEfinder interactively in bash
for f in /projects/cooper_research1/Wild_Sorghum_WGS/bam_wild/G*.bam; do
    name=$(basename $f| cut -f1 -d'.')
    mkdir $name
    cd $name

    # Cluster submission using Slurm
    sbatch -t '168:00:00' -N 1 --mem=48gb --ntasks-per-node=32 -o $name'_TIP'.%j --wrap="bash /projects/cooper_research1/TIP_Analysis/bin/TEfinder/TEfinder -alignment $f -fa $SORGUM_REF -gtf $TE_REF -te $TE_LIST"

## Single run - Not advisable for efficiency
# bash ~/bin/TEfinder/TEfinder -alignment $f -fa /scratch/sburkes/Sorghum/Sbicolor_454_v3.0.1.fa -gtf /scratch/sburkes/Results/results_Sbicolor_454_v3/repeatmasker_out/Sbicolor_454_v3.0.1.fa.out.gff -te ~/bin/TEfinder/TE_list.txt

    cd ..
done
```



## Analysis & Visualizations
### Value Counts of Insertions

To derive characteristics about detected  insertion sites, we're using a mix of bash and python with the pandas module (Last update:  8/31 7:38pm):


```R
# Aggregate insertion site bed files
cat *sites.bed > /nobackup/cooper_research/Shel/TIP_Analysis/merged_insertion-sites.bed

# Import libraries and bed file.
import pandas as pd
sb = pd.read_csv('Sbicolor_454_v3.0.1.bed', delimiter='\t', index_col=False, names=['Chromosome', 'Start','Stop', 'ID',5,6,7,8,9,10,'Superfamily',12,13,14])

# Get value_counts of unique TE superfamilies
sb['Superfamily'].value_counts()


# Formatted .bed for density visualization

```



### Visualization: Creating Plots using R [(Source)](https://www.biostars.org/p/69748/)

```R
# Check if ggplot2 is installed, if so, load it. If not, install and load it

if("ggplot2" %in% rownames(installed.packages())){
    library(ggplot2)
} else {
    install.packages("ggplot2")
    library(ggplot2)
}
    
# Import a text file with gene positions
# Column Headers: Chr, Position (No end or Gene name required)

genes <- read.table("genes.txt",sep="\t",header=T)
    
# make sure the chromosomes are ordered in the way you want them to appear in the plot

genes$chr <- with(genes, factor(chr, levels=paste("chr",c(1:22,"X","Y"),sep=""), ordered=TRUE))
    
# make a density plot of genes over the provided chromosomes (or scaffolds ...)

plottedGenes <- ggplot(genes) + geom_histogram(aes(x=pos),binwidth=1000000) + facet_wrap(~chr,ncol=2) + ggtitle("RefSeq genes density over human genome 19") + xlab("Genomic position (bins 1 Mb)") + ylab("Number of genes")
    
# Save it to an image
png("genes.png",width=1000,height=1500)
print(plottedGenes)
dev.off()
```

**Calculate Gene Density Per Kb And Plot Density Over Position For All Scaffolds Of A Draft Genome Using R**