## SNP calling on wild house mice populations

> Pull down data from [Harr et al.](https://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/)
```bash
## Note: bam files for North American (NH/VT) and South American (MAN) were kindly provided by members of the Nachman Lab

# M.m.domesticus
wget http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_domesticus/genomes_bam/*.bam
# M.m.musculus
wget http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_musculus/genomes_bam/*.bam
# M.m.castaneus
wget http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_castaneus/genomes_bam/*.bam

samtools sort -o ${Sample}_merge.sort.bam ${line}.bam
```

> Prepare for SNP calling with Picard
```bash
picard MarkDuplicates INPUT=${Sample}_merge.sort.bam OUTPUT=${Sample}_markdups.bam METRICS_FILE=${Sample}_metrics.txt

picard AddOrReplaceReadGroups I=${Sample}_markdups.bam O=${Sample}_markdups.rehead.bam RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM={Sample}

picard BuildBamIndex INPUT=${Sample}_markdups.rehead.bam

## HaplotypeCaller
gatk HaplotypeCaller -R Mus_musculus.GRCm38.dna.toplevel.fa -I ${Sample}_markdups.rehead.split.bam -ERC GVCF -stand-call-conf 20 -O ${Sample}_rawvariants.g.vcf.gz
```

> Filter SNPs
```bash
### Note: replace names of input files (ind1_rawvariants.g.vcf.gz) below with output of previous step HaplotypeCaller (e.g., *_rawvariants.g.vcf.gz)

# Combine files
gatk CombineGVCFs -R Mus_musculus.GRCm38.dna.toplevel.fa --variant ind1_rawvariants.g.vcf.gz --variant ind2_rawvariants.g.vcf.gz -O Combined_Ind.g.vcf.gz

## GenotypeGVCFs 
gatk --java-options \"-Xmx4g\" GenotypeGVCFs -R Mus_musculus.GRCm38.dna.toplevel.fa -V Combined_Ind.g.vcf.gz -O Combined_Ind.vcf.gz

# Select only SNPS and filter for low quality variants
gatk SelectVariants --reference Mus_musculus.GRCm38.dna.toplevel.fa --variant Combined_Ind.vcf.gz --select-type-to-include SNP --output Combined_Ind.SNPs.vcf.gz
gatk VariantFiltration --reference Mus_musculus.GRCm38.dna.toplevel.fa --variant Combined_Ind.SNPs.vcf.gz --filter-expression \"QD < 2.0 || FS > 60\" --filter-name \"SNPFilter\" --output Combined_Ind.SNPsfilt.vcf.gz
```
