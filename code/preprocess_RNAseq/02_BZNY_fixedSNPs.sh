#!/bin/bash

##Prepare script for genotyping multiple samples with GATK
##Requires GATK, samttols, picard in path

## Samples_for_genotyping.txt contains list of file prefixes
cat Samples_for_mapping_genotyping.txt | while read line
do
echo "
## Mapping to mm10 reference genome
bowtie2 -p 12 --very-sensitive -x MusGRCm38/GRCm38_68 -1 ${line}_L001_R1_001.fastq.gz,${line}_L002_R1_001.fastq.gz,${line}_L003_R1_001.fastq.gz,${line}_L004_R1_002.fastq.gz,${line}_L005_R1_002.fastq.gz -2 ${line}_L001_R2_001.fastq.gz,${line}_L002_R2_001.fastq.gz,${line}_L003_R2_001.fastq.gz,${line}_L004_R2_001.fastq.gz,${line}_L005_R2_001.fastq.gz -S ${line}.sam

samtools view -S -b ${line}.sam > ${line}.bam

samtools sort -o ${line}_merge.sort.bam ${line}.bam

## Prepare for SNP calling with picard
picard MarkDuplicates INPUT=${line}_merge.sort.bam OUTPUT=${line}_markdups.bam METRICS_FILE=${line}_metrics.txt

picard AddOrReplaceReadGroups I=${line}_markdups.bam O=${line}_markdups.rehead.bam RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM={Sample}

picard BuildBamIndex INPUT=${line}_markdups.rehead.bam

## HaplotypeCaller
gatk HaplotypeCaller -R Mus_musculus.GRCm38.dna.toplevel.fa -I ${line}_markdups.rehead.split.bam -ERC GVCF -stand-call-conf 20 -O ${line}_rawvariants.g.vcf.gz
"
done

###For next portion, replace names of input files (ind1_rawvariants.g.vcf.gz) below with output of previous step HaplotypeCaller (e.g., *_rawvariants.g.vcf.gz)
echo "
# Combine files
gatk CombineGVCFs -R Mus_musculus.GRCm38.dna.toplevel.fa --variant ind1_rawvariants.g.vcf.gz --variant ind2_rawvariants.g.vcf.gz -O Combined_Ind.g.vcf.gz

## GenotypeGVCFs 
gatk --java-options \"-Xmx4g\" GenotypeGVCFs -R Mus_musculus.GRCm38.dna.toplevel.fa -V Combined_Ind.g.vcf.gz -O Combined_Ind.vcf.gz

# Select only SNPS and filter for low quality variants
gatk SelectVariants --reference Mus_musculus.GRCm38.dna.toplevel.fa --variant Combined_Ind.vcf.gz --select-type-to-include SNP --output Combined_Ind.SNPs.vcf.gz
gatk VariantFiltration --reference Mus_musculus.GRCm38.dna.toplevel.fa --variant Combined_Ind.SNPs.vcf.gz --filter-expression \"QD < 2.0 || QUAL < 30.0 || FS > 200 || ReadPosRankSum < -20.0\" --filter-name \"SNPFilter\" --output Combined_Ind.SNPsfilt.vcf.gz

"
