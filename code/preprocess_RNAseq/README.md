### Trim and clean raw reads with [FastP v.0.19.6](https://github.com/OpenGene/fastp) [(Chen et al., 2018)](https://doi.org/10.1093/bioinformatics/bty560) ([runFastP.MAB.May21.sh]())
```
fastp \
        -i $seq\_R1.fastq \ # input read1 file
        -I $seq\_R2.fastq \ # input read2 file
        -o $seq\_R1_cleaned.fq \ # output read1 file
        -O $seq\_R2_cleaned.fq \ # output read 2 file
        -n 5 \ # read pair is discarded if number of N bases is >5
        -q 15 \ # minimum base quality score to keep
        -u 30 \ # percent of bases allowed to be less than q in a read
        --detect_adapter_for_pe \ # enable PE adapter trimming
        --cut_window_size=4 \
        --cut_mean_quality=15 \ # mean base score across the window required, or else trim the last base
        --length_required=25 \ # minimum read length to keep after trimming
        -j $seq\_report.json \ # output file name, JSON format
        -h $seq\_report.html \ # output file name, HTML format
        -w 2 # number of threads to use
```

### Align sequences to mouse genome via [STAR v.2.7.7a](https://github.com/alexdobin/STAR) ([IndexGenome.MAB.Jan21.sh]())
```
STAR --runThreadN 16 \
--runMode genomeGenerate \
--limitGenomeGenerateRAM 33524399488 \
--genomeDir /genomedirectory/ \
--genomeFastaFiles Mus_musculus.GRCm38.dna.toplevel.fa \
--sjdbGTFfile Mus_musculus.GRCm38.98.gtf \
--sjdbOverhang 149 # Readlength-1
```

> Align sequences with [STAR](https://github.com/alexdobin/STAR) ([STARalign.MAB.May21.sh]())
```
STAR --runMode alignReads \
--runThreadN 16 \
--genomeDir /genomedirectory \
--readFilesIn ${Sample}\_R1_cleaned.fq ${Sample}\_R2_cleaned.fq \
--outSAMtype BAM SortedByCoordinate \
--outFilterMultimapNmax 1 \ # max number of multiple alignments allowed for a read
--outFilterMismatchNmax 3 \ # max number of mismatches per pair
--outFileNamePrefix /mapped_reads/${Sample}
```

### Count reads with [HTseq v.0.11.0](https://htseq.readthedocs.io/en/release_0.11.1/index.html) ([HTSEQ_count.MAB.Jun21.sh]())
```
htseq-count -f bam --order=pos --stranded=no /mapped_reads/${Sample}_L001Aligned.sortedByCoord.out.bam Mus_musculus.GRCm38.98.gtf > /countfiles/${Sample}.L001.count

htseq-count -f bam --order=pos --stranded=no /mapped_reads/${Sample}_L004Aligned.sortedByCoord.out.bam Mus_musculus.GRCm38.98.gtf > /countfiles/${Sample}.L004.count

# merge lanes together (can also do this step at right at the beginning)
paste /countfiles/${Sample}.L001.count /countfiles/${Sample}.L004.count | awk -F' ' '{print $1"\t"$2+$4}' > /countfiles/${Sample}.count.merge
```

### Merge count files with [merge_tables.py](https://github.com/aiminy/SCCC-bioinformatics/blob/master/merge_tables.py) ([MergeCountTables.MAB.Jun21.sh])
```
python merge_tables.py \

/countfiles/all_compare.txt > /countfiles/all_compare_merged_counts.txt

# format of 'all_compare.txt':
# {Sample}.count.merge Pop_Trt_Sex_{Sample}  (e.g.: 002.count.merge BZrtM_002)

```





