#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=reorderBAM_slurm_job
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00

cd /projects/cooper_research1/TERRA_BAP/
mkdir -p /projects/cooper_research1/TERRA_BAP/BAM_by_ID/

## Get IDs of BAMs
#ls /projects/cooper_research1/TERRA_BAP/BAM_by_CHR/Chr01/*.bam | cut -d'/' -f7 | cut -d'.' -f1 > bamIDs.txt

## Build BAM list to loop through
bamARRAY=()
for f in /projects/cooper_research1/TERRA_BAP/BAM_by_CHR/Chr*; do
     bamARRAY+=($f)
done

## Loop through array to copy files
while read p; do
    for ELEMENT in ${bamARRAY[@]}; do
        mkdir -p /projects/cooper_research1/TERRA_BAP/BAM_by_ID/${p}
        cp -R -u -p ${ELEMENT}/${p}*.bam /projects/cooper_research1/TERRA_BAP/BAM_by_ID/${p}
    done
done </projects/cooper_research1/TERRA_BAP/bamIDs.txt
