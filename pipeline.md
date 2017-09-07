
# Overall steps:

The pipeline for a single paired-end lib contains 15 steps as follow:

1. report read length 
2. alignment 
3. flagstat bam 
4. post align: dedup bam 
5. post align: picard markedup 
6. post align: dedup bam (again) - final bam file 
7. post align: name sort bam 
8. post align: bam to bedpe 
9. bedpe to tagalign 
10. shift tagalign 
11. xcor subset sample
12. xcor calculation use subset 
13. macs2 peak calling 
14. filter peaks 
15. ataqc 

# Individual steps

## Alignment

```shell
bowtie2  -X2000 --mm --local | samtools view -Su /dev/stdin | samtools sort & index > xxx.PE2SE.bam &.bai
```

For bowtie2:

* Use memory-mapped I/O to load the index (--mm); 
* '-X2000' means maximum fragment length for valid paired-end alignments is 2000bp; 
* --local: a preset options mode, default as --sensitive-local, 

For samtools view:

* -S: ignore for compatibility with previous samtools versions 
* -u: uncompressed BAM outputs 

## Filter & deduplicate bam

``` shell
samtools view -F 1804 -f 2 -u -q 30 xxx.PE2SE.bam | sambamba sort -n  /dev/stdin -o /output_dir/xxx.PE2SE.dupmark.bam
```

1. Remove improper mapping marker (1804) & poor mapping score (<30) & output
 [u]ncompressed bam & [f] output fwd and rev. both mapped pairs 
2. Sort the bam by name (-n) and prepair for the deduplicating step 

``` shell
samtools fixmate -r xxx.PE2SE.dupmark.bam (tmp)  xxx.PE2SE.dupmark.bam.fixmate.bam (tmp) 
```

Fill in mate coordinate. ISIZE (insert size) and mate related flags from the name-sorted bam and remove secondary and ummapped reads (-r)

``` shell
 samtools view -F 1804 -f 2 -u xxx.PE2SE.dupmark.bam.fixmate.bam | sambamba sort  /dev/stdin -o xxx.PE2SE.filt.bam 
```


## Call peaks

``` shell
macs2 callpeak -t xxx.PE2SE.nodup.tn5.tagAlign.gz -f BED \
-n xxx.PE2SE.nodup.tn5.pf" -g "hs" -p 0.01 --nomodel \ 
--shift -75 --extsize 150 -B --SPMR --keep-dup all --call-summits 
```

Sort by Col8 in descending order and replace long peak names in Column 4 with Peak_

``` shell
sort -k 8gr,8gr xxx.PE2SE.nodup.tn5.pf"_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | gzip -nc > xxx.PE2SE.nodup.tn5.pf.narrowPeak.gz
```

``` shell
 macs2 bdgcmp -t xxx.PE2SE.nodup.tn5.pf"_treat_pileup.bdg -c xxx.PE2SE.nodup.tn5.pf"_control_lambda.bdg \
 --o-prefix xxx.PE2SE.nodup.tn5.pf" -m FE
 
slopBed -i xxx.PE2SE.nodup.tn5.pf"_FE.bdg -g hg38.chrom.sizes -b 0 | bedClip stdin hg38.chrom.sizes xxx.PE2SE.nodup.tn5.pf.fc.signal.bedgraph

sort -k1,1 -k2,2n xxx.PE2SE.nodup.tn5.pf.fc.signal.bedgraph > xxx.PE2SE.nodup.tn5.pf.fc.signal.srt.bedgraph


bedGraphToBigWig xxx.pf.fc.signal.srt.bedgraph hg38.chrom.sizes xxx.PE2SE.nodup.tn5.pf.fc.signal.bigwig

```

## Lib QCs

Some concerpts:

* insert size 
* fragment distribution 

### tss enrichment caclculation

* Calculated by using the final bam file 
* Extended TSS to -/+2kb 
* Use metaseq package to create [BamSignal class](https://pythonhosted.org/metaseq/autodocs/metaseq._genomic_signal.BamSignal.html), and caclulated coverageover TSS features which stores in a  (length(features)*bins) NumPy array
* Shifted the bam file to half of the read length in the 5' direction 
* Reversed the promoters on the minus strand 
* Use normalization method from Greenleaf et al. 2013:
  * background average noise is to use averaged coverage of 100bps at both ends 
  * enrichment = coverage / background average noise



# Reference
1. [ENCODE Standards ](https://www.encodeproject.org/atac-seq/)

  


