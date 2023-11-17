module load cellranger/7.2.0

#Define file path where raw FASTQ files from sequencing are saved
FASTQ_PATH=


##Download human reference genome and annotation from Ensembl
wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz

##Uncompress files
zcat Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > GRCh38.110.fa
zcat Homo_sapiens.GRCh38.110.gtf.gz > GRCh38.110.gtf

##Concatenate the human and IAV Cal09 reference genomes
cp GRCh38.110.fa GRCh38.110_cal09.fa
cat cal09.fasta >> GRCh38.110_cal09.fa

##Filter the human gtf file
cellranger mkgtf \
  GRCh38.110.gtf \
  GRCh38.110.filtered.gtf \
  --attribute=gene_biotype:protein_coding

##Add the Cal09 GTF to the end of the filtered GRCh38 gtf file
cp GRCh38.110.filtered.gtf GRCh38.110.filtered_cal09.gtf
cat cal09.gtf >> GRCh38.110.filtered_cal09.gtf


##Create cell ranger genome reference
cellranger mkref \
--genome=ref \
--fasta=GRCh38.110_cal09.fa \
--genes=GRCh38.110.filtered_cal09.gtf

#Cleanup files, re-organize output directory
mv ref/* ./
rm -rf ref/

#Cellranger count the Cal09 infected sample
cellranger count --id=cal_out \
                 --transcriptome=./ref \
                 --fastqs= $FASTQ_PATH \
                 --sample=Langlois_016_S2_FL
                 
#Cellranger count the mock infected sample
cellranger count --id=mock_out \
                 --transcriptome=./ref \
                 --fastqs= $FASTQ_PATH \
                 --sample=Langlois_016_S5_FL