#!/bin/bash
##Obtain count data for parental (NY, BZ) RNA-seq reads
##Requires FastP, STAR, HTseq-count

cat filenames.txt | while read -r line
do
echo "

##trim and clean RNAseq Reads
fastp -i ${line}\_R1.fastq -I ${line}\_R2.fastq -o ${line}\_R1_cleaned.fq -O ${line}\_R2_cleaned.fq \
      -n 5 -q 15 -u 30 --detect_adapter_for_pe --cut_window_size=4 --cut_mean_quality=15 --length_required=25 \
      -j ${line}\_report.json -h ${line}\_report.html -w 2

##index genome
STAR --runThreadN 16 --runMode genomeGenerate --limitGenomeGenerateRAM 33524399488 --genomeDir genome_index_STAR_mm10 \
	 --genomeFastaFiles Mus_musculus.GRCm38.dna.toplevel.fa --sjdbGTFfile Mus_musculus.GRCm38.98.gtf --sjdbOverhang 149

##align RNAseq reads
STAR --runMode alignReads --runThreadN 16 --genomeDir genome_index_STAR_mm10 --readFilesIn ${line}\_R1_cleaned.fq ${line}\_R2_cleaned.fq \
	 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 --outFilterMismatchNmax 3 --outFileNamePrefix ${line}

##produce count files
python -m HTSeq.scripts.count -f bam --order=pos --stranded=no ${line}_L001Aligned.sortedByCoord.out.bam Mus_musculus.GRCm38.98.gtf > ${line}.L001.count

python -m HTSeq.scripts.count -f bam --order=pos --stranded=no ${line}_L004Aligned.sortedByCoord.out.bam Mus_musculus.GRCm38.98.gtf > ${line}.L004.count

##merge lanes together
paste ${line}.L001.count ${line}.L004.count | awk -F' ' '{print $1"\t"$2+$4}' > ${line}.count.merge

##merge count files together
python merge_tables.py all_parents.txt > all_parents_counts.txt

"
done
