### Four main analyses are outlined here:
1) Counting reads for parental samples
2) Identifying fixed SNPs between Brazil and New York
3) Counting alleleic reads in F1 hybrids
4) Identifying *cis* and *trans* regulatory divergence

### 1) Parental Read Counts

> Trim and clean raw reads with [FastP v.0.19.6](https://github.com/OpenGene/fastp)
```bash
fastp -i ${Sample}\_R1.fastq -I ${Sample}\_R2.fastq -o ${Sample}\_R1_cleaned.fq -O ${Sample}\_R2_cleaned.fq \
      -n 5 -q 15 -u 30 --detect_adapter_for_pe --cut_window_size=4 --cut_mean_quality=15 --length_required=25 \
      -j ${Sample}\_report.json -h ${Sample}\_report.html -w 2
```

> Index mouse genome via [STAR v.2.7.7a](https://github.com/alexdobin/STAR)
```bash
STAR --runThreadN 16 --runMode genomeGenerate --limitGenomeGenerateRAM 33524399488 --genomeDir genome_index_STAR_mm10 \
     --genomeFastaFiles Mus_musculus.GRCm38.dna.toplevel.fa --sjdbGTFfile Mus_musculus.GRCm38.98.gtf --sjdbOverhang 149
```

> Align sequences with [STAR v.2.7.7a](https://github.com/alexdobin/STAR)
```bash
STAR --runMode alignReads --runThreadN 16 --genomeDir genome_index_STAR_mm10 \
     --readFilesIn ${Sample}\_R1_cleaned.fq ${Sample}\_R2_cleaned.fq --outSAMtype BAM SortedByCoordinate \
     --outFilterMultimapNmax 1 --outFilterMismatchNmax 3 --outFileNamePrefix ${Sample}
```

> Count reads with [HTseq v.0.11.0](https://htseq.readthedocs.io/en/release_0.11.1/index.html)
```python
python -m HTSeq.scripts.count -f bam --order=pos --stranded=no ${Sample}_L001Aligned.sortedByCoord.out.bam Mus_musculus.GRCm38.98.gtf > ${Sample}.L001.count

python -m HTSeq.scripts.count -f bam --order=pos --stranded=no ${Sample}_L004Aligned.sortedByCoord.out.bam Mus_musculus.GRCm38.98.gtf > ${Sample}.L004.count

# merge lanes together (can also do this step at right at the beginning)
paste ${Sample}.L001.count ${Sample}.L004.count | awk -F' ' '{print $1"\t"$2+$4}' > ${Sample}.count.merge
```

> Merge count files with [merge_tables.py](https://github.com/aiminy/SCCC-bioinformatics/blob/master/merge_tables.py)
```python
python merge_tables.py all_parents.txt > all_parents_counts.txt

# format of 'all_parents.txt':
# {Sample}.count.merge Pop_Trt_Sex_{Sample}  (e.g.: 002.count.merge BZrtM_002)
```
###### Note: _all_parents_counts.txt_ is provided in data/raw/ReadCounts/all_parents_counts.txt ######


### 2) Fixed SNPs between BZ and NY

> Map genomic reads to mouse reference genome via [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
```bash


samtools sort -o ${Sample}_merge.sort.bam ${Sample}_merge.bam
````

> Mark duplicates with [Picard](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-)
```java
picard MarkDuplicates INPUT=${Sample}_merge.sort.bam OUTPUT=${Sample}_markdups.bam METRICS_FILE=${Sample}_metrics.txt

picard BuildBamIndex INPUT=${Sample}_markdups.bam

picard AddOrReplaceReadGroups I=${Sample}_markdups.bam O=${Sample}_markdups.rehead.bam \
       RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${Sample}
```

> Joint genotyping via GATK [HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) and [GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/360037057852-GenotypeGVCFs)
```bash
gatk HaplotypeCaller -R Mus_musculus.GRCm38.dna.toplevel.fa -I ${Sample}_markdups.rehead.split.bam -ERC GVCF \
     -stand-call-conf 20 -O ${Sample}_rawvariants.g.vcf.gz

#Combine files
gatk CombineGVCFs -R Mus_musculus.GRCm38.dna.toplevel.fa --variant MANA_rawvariants.g.vcf.gz \
     --variant SARA_rawvariants.g.vcf.gz -O Combined_BZ_NY.g.vcf.gz

gatk --java-options \"-Xmx4g\" GenotypeGVCFs -R Mus_musculus.GRCm38.dna.toplevel.fa -V Combined_BZ_NY.g.vcf.gz \
     -O Combined_BZ_NY.vcf.gz
```

> Select and filter variants via GATK [SelectVariants](https://gatk.broadinstitute.org/hc/en-us/articles/360037055952-SelectVariants) and [VariantFiltration](https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration)
```bash
gatk SelectVariants --reference Mus_musculus.GRCm38.dna.toplevel.fa --variant Combined_BZ_NY.vcf.gz \
     --select-type-to-include SNP --output Combined_BZ_NY.SNPs.vcf.gz

gatk VariantFiltration --reference Mus_musculus.GRCm38.dna.toplevel.fa --variant Combined_BZ_NY.SNPs.vcf.g vcf gz \
     --filter-expression \"QD < 2.0 || QUAL < 30.0 || FS > 200.0 || ReadPosRankSum < -20.0\" --filter-name \"SNPFilter\" \
     --output Combined_BZ_NY.SNPsfilt.vcf.gz
```

### 3) Allelic Read Counts

###### Note: F1 sequences were already trimmed, cleaned, and mapped as described above under '1) Parental Read Counts'. Here, we are remapping F1 sequences:

> Align F1 sequences with [STAR v.2.7.7a](https://github.com/alexdobin/STAR)
```bash
STAR --runMode alignReads --runThreadN 16 --genomeDir genome_index_STAR_mm10 \
     --readFilesIn ${Sample}_R1.trim.fastq ${Sample}_R2.trim.fastq --outSAMtype BAM SortedByCoordinate \
     --waspOutputMode SAMtag --outSAMattributes NH HI AS nM NM MD vA vG vW \
     --varVCFfile BZ_NY_filteredhetcalls.vcf --outFilterMultimapNmax 1 \
     --outFilterMismatchNmax 3 --outFileNamePrefix ${Sample}.STAR_ASE

samtools view -h -o ${Sample}.STAR_ASEAligned.sortedByCoord.out.sam ${Sample}.STAR_ASEAligned.sortedByCoord.out.bam 
```
###### Note: _BZ_NY_filteredhetcalls.vcf_ was produced in step '(2) Identifying fixed SNPs between BZ and NY'

> Filter files based on [WASP]() flags
```bash
# Sort reads into allele-specific pools (NY, BZ)

#vW:i:1 means alignment passed WASP filtering, and all other values mean it did not pass:
grep \"vW:i:1\" ${Sample}.STAR_ASEAligned.sortedByCoord.out.sam > ${Sample}.STAR_ASEAligned.sortedByCoord.out.filter.sam
# Note: geno1 and geno2 refer to NY and BZ reads
grep \"vW:i:1\" ${Sample}.STAR_ASEAligned.sortedByCoord.out.sam | grep \"vA:B:c,1\" \
     | awk -F' ' '{if(\$18 ~ /2/){}else{print}}' > ${Sample}.geno1.sam
grep \"vW:i:1\" ${Sample}.STAR_ASEAligned.sortedByCoord.out.sam | grep \"vA:B:c,2\" \
     | awk -F' ' '{if(\$18 ~ /1/){}else{print}}' > ${Sample}.geno2.sam

grep \"@\" ${Sample}.geno2.sam > ${Sample}.head.sam

cat ${Sample}.head.sam ${Sample}.geno1.sam > ${Sample}.geno1.withhead.sam
cat ${Sample}.head.sam ${Sample}.geno2.sam > ${Sample}.geno2.withhead.sam
cat ${Sample}.head.sam ${Sample}.STAR_ASEAligned.sortedByCoord.out.filter.sam > ${Sample}.bothwhead.sam

samtools view -S -b ${Sample}.geno1.withhead.sam > ${Sample}.geno1.withhead.bam
samtools view -S -b ${Sample}.geno2.withhead.sam > ${Sample}.geno2.withhead.bam
samtools view -S -b ${Sample}.bothwhead.sam > ${Sample}.bothwhead.bam

samtools sort ${Sample}.geno1.withhead.bam -o ${Sample}.geno1.withhead.sorted.bam
samtools sort ${Sample}.geno2.withhead.bam -o ${Sample}.geno2.withhead.sorted.bam
samtools sort ${Sample}.bothwhead.bam -o ${Sample}.bothwhead.sorted.bam
````

> Produce count files via HTseq
```python
python -m HTSeq.scripts.count -f bam -s no \
       -r pos ${Sample}.geno1.withhead.sorted.bam Mus_musculus.GRCm38.98.gtf >  ${Sample}.geno1.withhead.sorted.counts

python -m HTSeq.scripts.count -f bam -s no \
       -r pos ${Sample}.geno2.withhead.sorted.bam Mus_musculus.GRCm38.98.gtf >  ${Sample}.geno2.withhead.sorted.counts
```

> Assign reads to specific genotypes via [GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-)
```java
java -jar picard-2.20.7/picard.jar AddOrReplaceReadGroups \
     I=${Sample}.geno1.withhead.sorted.bam O=${Sample}.geno1.withhead.plusreads.sorted.withgroup.bam \
     RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${Sample}_geno1

java -jar picard-2.20.7/picard.jar AddOrReplaceReadGroups \
     I=${Sample}.geno2.withhead.sorted.bam O=${Sample}.geno2.withhead.plusreads.sorted.withgroup.bam \
     RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${Sample}_geno2

java -jar picard-2.20.7/picard.jar AddOrReplaceReadGroups I=${Sample}.bothwhead.sorted.bam \
     O=${Sample}.bothwhead.withgroup.bam RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${Sample}_both
```

> Produce count tables per SNP via [GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360037428291-ASEReadCounter)
```bash
gatk ASEReadCounter -I ${line}.bothwhead.sorted.bam O=${line}.bothwhead.withgroup.bam \
     -R Mus_musculus.GRCm38.dna.toplevel.fa -V BR_NY_hetcalls.vcf.gz \
     -O BZ_NY.bothgeno.table

gatk ASEReadCounter -I ${line}.geno1.withhead.plusreads.sorted.withgroup.bam \
     -R Mus_musculus.GRCm38.dna.toplevel.fa -V BR_NY_hetcalls.vcf.gz \
     -O ${line}.geno1.table

gatk ASEReadCounter -I ${line}.geno2.withhead.plusreads.sorted.withgroup.bam \
     -R Mus_musculus.GRCm38.dna.toplevel.fa -V BR_NY_hetcalls.vcf.gz \
     -O ${line}.geno2.table

gatk ASEReadCounter -I ${line}.bothwhead.sorted.bam -R Mus_musculus.GRCm38.dna.toplevel.fa \
     -V BR_NY_hetcalls.vcf.gz -O ${line}.both_alleles.table
```

### 4) Patterns of *cis* and *trans*

> Read in count tables
```R
counts <- read.csv("counts_BAT_males.txt",sep="\t",header=TRUE,row.names=1)
condTable <- read.csv("condTable_BAT_males.txt",sep="\t",header=TRUE)
counts2<-as.matrix(counts)
```

> Identify *trans* effects using LRT (contrasting log2 Fold Change between F1's and parents)
```R
#v1
design = ~F1_Parent + population + temp +population:F1_Parent
dds <- DESeqDataSetFromMatrix(counts2, condTable, design)
dds <- DESeq(dds, test="LRT", reduced= ~ population + F1_Parent + temp)

#v2
design = ~F1_Parent + population + population:F1_Parent
dds <- DESeqDataSetFromMatrix(counts2, condTable, design)
dds <- DESeq(dds, test="LRT", reduced= ~ population + F1_Parent)
res <- results(dds)
```

> Identify *cis* effects
```R
# make sure there is a pseudo-sample column where samples are rep'd across treatments of interest
# pseudo-sample is sample trick from Michael Love, so that temp info can be incorporated

dds <- DESeqDataSetFromMatrix(counts2, condTable, design)

#m = # of individuals
design = ~temp + temp:pseudo_sampl + temp:allele
dds <- DESeqDataSetFromMatrix(counts2, condTable, design)
m <- 6

#don't estimate size factors here b/c this is ASE
sizeFactors(dds) <- rep(1, 2*m)
dds <- DESeq(dds, fitType="parametric")
resultsNames(dds)
```

> Sort categories by comparing log2 fold change and padj values between parents, F1's, and LRT (for trans component)
```bash
# Note: we adjusted FDR for cases where genes were tested by DESeq2 for one test but not the other
# (i.e., FDR estimates for parents but not hybrid) -- this is done with p.adjust, method="BH"

#0.05 cutoff
awk -F' ' '{if($3<0.05 && $5<0.05 && $7>0.05){print $0,"CIS_only"}else{print}}' |
awk -F' ' '{if($3>0.05 && $5<0.05 && $7<0.05){print $0,"Compensatory"}else{print  $0}}'|
awk -F' ' '{if($3>0.05 && $5>0.05 && $7>0.05){print $0,"Conserved"}else{print $0}}' |
awk -F' ' '{if($3<0.05 && $5>0.05 && $7<0.05){print $0,"Trans"}else{print  $0}}' |
awk -F' ' '{if( ($3<0.05 && $5<0.05 && $7<0.05) && ( ($2>0 && $4>0) || ($2<0 && $4<0) ) && ( sqrt($2^2)>sqrt($4^2)) ){print $0,"cis+trans_same"}else{print}}' |
awk -F' ' '{if( ($3<0.05 && $5<0.05 && $7<0.05) && ( ($2>0 && $4>0) || ($2<0 && $4<0) ) && ( sqrt($2^2)<sqrt($4^2)) ){print $0,"cis+trans_opp"}else{print}}' |
awk -F' ' '{if( ($3<0.05 && $5<0.05 && $7<0.05) && ( ($2>0 && $4<0) || ($2<0 && $4>0) )){print $0,"cisxtrans"}else{print}}' |
awk -F' ' '{if($3=="NA" || $5=="NA" || $7=="NA"){print $1,$2,$3,$4,$5,$6,$7,"NoPower"}else{print}}' |
awk -F' ' '{if($8==""){print $0,"Ambiguous"}else{print}}' 

# add the following to condense cis+trans groups for plotting in R
 | awk -F' ' '{if($8=="Conserved" || $8=="Ambiguous"){print $1,$2,$3,$4,$5,$6,$7,"Ambiguous_orConserved"}else{print}}'  > categories.txt

```

