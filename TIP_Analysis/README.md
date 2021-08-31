# TIPs Detection Workflow in Sorghum

## Pre-reqs

- Rmsk2bed
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
# Change directory to working directory of TEFinder
cd /users/sburkes/bin/TEfinder
    
# Location of the reference genome fasta file
SORGUM_REF='/scratch/sburkes/Sorghum/Sbicolor_454_v3.0.1.fa'
    
# Removes simple repeats from RepeatMasker .out
grep -v -iE '(Motif\:[ATGC]+\-rich)|(Motif\:\([ATGC]+\)n)' /scratch/sburkes/Results/results_Sbicolor_454_v3/repeatmasker_out/Sbicolor_454_v3.0.1.fa.out > TEs.gtf
    
# Converts .out from RepeatMasker into .bed file.
rmsk2bed < TEs.gtf bedops --merge - > Sbicolor_454_v3.0.1.bed
    
## Commented out 08/12/21 - Limits to LTR/Copia elements, truncates to the first element entry in .bed
# cat Sbicolor_454_v3.0.1.bed | cut -d, -f11 | grep "LTR/Copia" > truncated.bed
    
cat Sbicolor_454_v3.0.1.bed > truncated.bed
    
# Removes asterisk characters
tr -d '*' < truncated.bed > Sbicolor_truncated_rmchar.bed
    
# Converts .bed coordinates from .bed to fasta
bedtools getfasta -fi $SORGUM_REF -name -bed Sbicolor_truncated_rmchar.bed > Sbicolor_454_v3.0.1_TE_truncated.fa
    
sed 's/::.*//' Sbicolor_454_v3.0.1_TE_truncated.fa > Sbicolor_454_v3.0.1_TE_truncated_format.fa
    
# Outputs fasta entries into .txt file
grep -e ">" Sbicolor_454_v3.0.1_TE_truncated_format.fa | awk 'sub(/^>/, "")' >> TE_list.txt
        
```

## Detect TE insertions

Once the setup and gathering of files is complete, we can proceed with the discovery of non-reference TE insertions using the following script:

```bash
# Change directory to working directory of TEFinder
cd /users/sburkes/bin/TEfinder

# Execute TEFinder script with prior generated files
bash /users/sburkes/bin/TEfinder/TEfinder -alignment $f -fa /scratch/sburkes/Sorghum/Sbicolor_454_v3.0.1.fa -gtf /scratch/sburkes/Results/results_Sbicolor_454_v3/repeatmasker_out/Sbicolor_454_v3.0.1.fa.out.gff -te /users/sburkes/bin/TEfinder/TE_list.txt

# For the sake of time, this can be looped like so:
# Running TEfinder interactively in bash
for f in /projects/cooper_research1/Wild_Sorghum_WGS/bam_wild/G*.bam; do
    name=$(basename $f| cut -f1 -d'.')
    printf $name
    mkdir $name
    cd $name

  # Cluster submission using Slurm
  sbatch -t '72:00:00' -N 1 --mem=48gb --ntasks-per-node=32 -o $name'_TIP'.%j --wrap="bash /users/sburkes/bin/TEfinder/TEfinder -alignment $f -fa /scratch/sburkes/Sorghum/Sbicolor_454_v3.0.1.fa -gtf /scratch/sburkes/Results/results_Sbicolor_454_v3/repeatmasker_out/Sbicolor_454_v3.0.1.fa.out.gff -te /users/sburkes/bin/TEfinder/TE_list.txt"
    cd ..
    done

```




## Analysis
To derive characteristics about detected  insertion sites, we're using a mix of bash and python with the pandas module (Last update:  8/31 7:38pm):


```R
# Aggregate insertion site bed files
cat *sites.bed > /nobackup/cooper_research/Shel/TIP_Analysis/merged_insertion-sites.bed

# Import libraries and bed file.
import pandas as pd
sb = pd.read_csv('Sbicolor_454_v3.0.1.bed', delimiter='\t', index_col=False, names=['Chromosome', 'Start','Stop', 'ID',5,6,7,8,9,10,'Superfamily',12,13,14])

# Get value_counts of unique TE superfamilies
sb['Superfamily'].value_counts()
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
# Column Headers: Chr, Position (no end or gene name required)

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