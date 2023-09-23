
################################################################################
# Construction of tree shrew reference genome
################################################################################

# STEP 1: download genomeannotation.gtf and genome.fa for tree shrew
# http://www.treeshrewdb.org/download.html

# STEP 2: MAKE GTF FILE
cellranger mkgtf \
  TS_3.0.genomeannotation.gtf \
  TS_3.0.genomeannotation.filtered.gtf \
  --attribute=gene_biotype:protein_coding

# STEP 3: MAKE GTF FILE
cellranger mkref \
  --genome=TS_3.0_genome \
  --fasta=TS_3.0_genome.fa \
  --genes=TS_3.0.genomeannotation.filtered.gtf

# STEP4: RUN Cellranger
nohup cellranger count --id=sampleID \
     --transcriptome=/home/fingerstyle/TS_3.0_genome \
     --fastqs=/home/fingerstyle/TS-1S \
     --sample=sampleID \
     --localcores=10 \
     --localmem=64  &
