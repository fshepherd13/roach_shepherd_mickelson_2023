#!/bin/bash

IN_DIR=$1 #Directory containing the segment-specific bam files, combined by replicate
OUT_FILE=$2 #Assign output file path provided by snakemake rule (path/to/sample_cmratio.txt)

cPB2=$(samtools view ${IN_DIR}/PB2_all.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAACG"|wc -l) 
cPB1=$(samtools view ${IN_DIR}/PB1_all.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAATG"|wc -l)
cPA=$(samtools view ${IN_DIR}/PA_all.bam| awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAAGT"|wc -l)
cNP=$(samtools view ${IN_DIR}/NP_all.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAATA"|wc -l)
cHA=$(samtools view ${IN_DIR}/HA_all.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAACA"|wc -l)
cNA=$(samtools view ${IN_DIR}/NA_all.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAACT"|wc -l)
cM=$(samtools view ${IN_DIR}/M_all.bam| awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAACT"|wc -l)
cNS=$(samtools view ${IN_DIR}/NS_all.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAACA"|wc -l)

mPB2=$(samtools view ${IN_DIR}/PB2_all.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAAA"|wc -l)
mPB1=$(samtools view ${IN_DIR}/PB1_all.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAAAA"|wc -l)
mPA=$(samtools view ${IN_DIR}/PA_all.bam| awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAAA"|wc -l)
mNP=$(samtools view ${IN_DIR}/NP_all.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAAA"|wc -l)
mHA=$(samtools view ${IN_DIR}/HA_all.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAAA"|wc -l)
mNA=$(samtools view ${IN_DIR}/NA_all.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAAA"|wc -l)
mM=$(samtools view ${IN_DIR}/M_all.bam| awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAAA"|wc -l)
mNS=$(samtools view ${IN_DIR}/NS_all.bam|awk 'BEGIN {FS="\t"}; {print $10}'|grep "AAAAAAA"|wc -l)

echo -e "gene\tmRNA\tcRNA" > $OUT_FILE
echo -e "PB2\t$mPB2\t$cPB2" >> $OUT_FILE
echo -e "PB1\t$mPB1\t$cPB1" >> $OUT_FILE
echo -e "PA\t$mPA\t$cPA" >> $OUT_FILE
echo -e "HA\t$mHA\t$cHA" >> $OUT_FILE
echo -e "NP\t$mNP\t$cNP" >> $OUT_FILE 
echo -e "NA\t$mNA\t$cNA" >> $OUT_FILE
echo -e "M\t$mM\t$cM" >> $OUT_FILE
echo -e "NS\t$mNS\t$cNS" >> $OUT_FILE

awk 'NR==1{$4="mRNA:total_pos_RNA"}NR>1{if($2+$3==0) $4 = "N/A"; else $4=$2/($2+$3)}{print $1 "\t" $2 "\t" $3 "\t" $4}' $OUT_FILE > ${IN_DIR}/testfile.tmp && mv ${IN_DIR}/testfile.tmp $OUT_FILE