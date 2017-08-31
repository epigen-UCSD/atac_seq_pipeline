# Alignment 
```shell
bowtie2  -X2000 --mm --local | samtools view -Su /dev/stdin | samtools sort & index 
```

# Filter bam 
```shell
samtools view -F 1804 -f 2 -u /home/zhc268/scratch/outputs/JYH_54/align/rep1/JYH_54_R1.fastq.bz2.PE2SE.dupmark.bam.fixmate.bam | sambamba sort -t 8 /dev/stdin -o /home/zhc268/scratch/outputs/JYH_54/align/rep1/JYH_54_R1.fastq.bz2.PE2SE.filt.bam
```


# Call peaks 
```shell
macs2 callpeak -t $tag -f BED -n "$prefix" -g "$gensz" -p $pval_thresh --nomodel --shift -$shiftsize --extsize $smooth_window -B --SPMR --keep-dup all --call-summits $extra_param_macs2
```
Sort by Col8 in descending order and replace long peak names in Column 4 with Peak_<peakRank>
```shell
sort -k 8gr,8gr "$prefix"_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | gzip -nc > $peakfile
```

