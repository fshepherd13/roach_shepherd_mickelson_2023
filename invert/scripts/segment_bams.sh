#!/bin/bash
#SBATCH -J "IAV segment bams"
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=30:00:00
#SBATCH -p amdsmall
#SBATCH --mem=5GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sheph085@umn.edu

#Script to extract the segment-specific full-length bam file from select nHTBE/A549 infection data. Bam files will be visualized in Geneious to check for coverage consistency
#Samples analyzed:
#JF291_5 = nHTBE, ciliated, 12hr, 2 MOI, replicate 1
#JF291_13 = nHTBE, ciliated, 6hr, 2 MOI, replicate 1
#JF291_25 = nHTBE, secretory, 12hr, 2 MOI, replicate 1
#JF291_33 = nHTBE, secretory, 6hr, 2 MOI, replicate 1
#JF291_45 = A549, 12hr, 2 MOI, replicate 1
#JF291_53 = A549, 6hr, 2 MOI, replicate 1


module load samtools

cd ../results/invert/
for F in "JF291_5"; do
    cd ${F}_*/bam
    samtools view ${F}*_fwd.bam Cal_PB2 -o PB2_fwd.bam
    samtools view ${F}*_rev.bam Cal_PB2 -o PB2_rev.bam
    samtools view ${F}*_fwd.bam Cal_PB1 -o PB1_fwd.bam
    samtools view ${F}*_rev.bam Cal_PB1 -o PB1_rev.bam
    samtools view ${F}*_fwd.bam Cal_PA -o PA_fwd.bam
    samtools view ${F}*_rev.bam Cal_PA -o PA_rev.bam
    samtools view ${F}*_fwd.bam Cal_HA -o HA_fwd.bam
    samtools view ${F}*_rev.bam Cal_HA -o HA_rev.bam
    samtools view ${F}*_fwd.bam Cal_NP -o NP_fwd.bam
    samtools view ${F}*_rev.bam Cal_NP -o NP_rev.bam
    samtools view ${F}*_fwd.bam Cal_NA -o NA_fwd.bam
    samtools view ${F}*_rev.bam Cal_NA -o NA_rev.bam    
    samtools view ${F}*_fwd.bam Cal_M -o M_fwd.bam
    samtools view ${F}*_rev.bam Cal_M -o M_rev.bam
    samtools view ${F}*_fwd.bam Cal_NS -o NS_fwd.bam
    samtools view ${F}*_rev.bam Cal_NS -o NS_rev.bam
    cd ../../
done