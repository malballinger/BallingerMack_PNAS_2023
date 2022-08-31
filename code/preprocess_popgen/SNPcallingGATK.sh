#!/bin/bash
##Prepare script for genotyping multiple samples with GATK
##Requires GATK, samttols, picard in path

cat Samples_for_genotyping.txt | while read line
do
echo "
samtools sort -o ${line}_merge.sort.bam ${line}_merge.bam

picard MarkDuplicates INPUT=${line}_merge.sort.bam OUTPUT=${line}_markdups.bam METRICS_FILE=${line}_metrics.txt

picard AddOrReplaceReadGroups I=${line}_markdups.bam O=${line}_markdups.rehead.bam RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${line}

picard BuildBamIndex INPUT=${line}_markdups.rehead.bam

gatk HaplotypeCaller -R Mus_musculus.GRCm38.dna.toplevel.fa -I ${line}_markdups.rehead.split.bam -ERC GVCF -stand-call-conf 20 -O ${line}_rawvariants.g.vcf.gz
"
done

echo "
#Combine files
gatk CombineGVCFs -R Mus_musculus.GRCm38.dna.toplevel.fa --variant indX_rawvariants.g.vcf.gz --variant indY_rawvariants.g.vcf.gz -O Combined_Ind.g.vcf.gz

gatk --java-options \"-Xmx4g\" GenotypeGVCFs -R Mus_musculus.GRCm38.dna.toplevel.fa -V Combined_Ind.g.vcf.gz -O Combined_Ind.vcf.gz

#Select only SNPS and filter for low quality variants
gatk SelectVariants --reference Mus_musculus.GRCm38.dna.toplevel.fa --variant Combined_Ind.vcf.gz --select-type-to-include SNP --output Combined_Ind.SNPs.vcf.gz
gatk VariantFiltration --reference Mus_musculus.GRCm38.dna.toplevel.fa --variant Combined_Ind.SNPs.vcf.gz --filter-expression \"QD < 2.0 || FS > 60.0\" --filter-name \"SNPFilter\" --output Combined_Ind.SNPsfilt.vcf.gz
"
done
