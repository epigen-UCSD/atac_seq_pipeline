ATAC-Seq  Pipeline in Epigen-UCSD
===================================================

This pipeline is largely based on [Encode pipeline](./README_legend.md), with few small changes.: 

1. Enabled handling `fastq.bz` files 
2. Enabled to run by lib 
3. Configged to the epigen-group@TSCC 


More detailed documentation can be found on [this gitbook](https://biomystery.gitbooks.io/ataq_seq/). 

# List of all parameters


```
== atac pipeline settings
	-type <string>                   : Type of the pipeline. atac-seq or chip-seq (default: atac-seq).
	-chip_seq <bool>                 : Chip-Seq (no tn5 shifting).
	-trimmed_fastq <bool>            : Skip fastq-trimming stage.
	-align <bool>                    : Align only (no MACS2 peak calling or IDR or ataqc analysis).
	-subsample_xcor <string>         : # reads to subsample for cross corr. analysis (default: 25M).
	-subsample <string>              : # reads to subsample exp. replicates. Subsampled tagalign will be used for steps downstream (default: 0; no subsampling).
	-true_rep <bool>                 : No pseudo-replicates.
	-no_ataqc <bool>                 : No ATAQC
	-no_xcor <bool>                  : No Cross-correlation analysis.
	-csem <bool>                     : Use CSEM for alignment.
	-smooth_win <string>             : Smoothing window size for MACS2 peak calling (default: 150).
	-idr_thresh <real>               : IDR threshold : -log_10(score) (default: 0.1).
	-old_trimmer <bool>              : Use legacy trim adapters (trim_galore and trimAdapter.py).
	-ENCODE3 <bool>                  : Force to use parameter set (-smooth_win 73 -idr_thresh 0.05 -multimapping 4) for ENCODE3.
	-ENCODE <bool>                   : Force to use parameter set (-smooth_win 73 -idr_thresh 0.05 -multimapping 4) for ENCODE.
	-no_browser_tracks <bool>        : Disable generation of genome browser tracks (workaround for bzip2 shared library issue).
	-overlap_pval_thresh <real>      : p-val threshold for overlapped peaks (default: 0.01).
	-macs2_pval_thresh <real>        : MACS2 p-val threshold for calling peaks (default: 0.1).
	-macs2_pval_thresh_bw <real>     : MACS2 p-val threshold for generating BIGWIG signal tracks (default: 0.1).
	-enable_idr <bool>               : Enable IDR on called peaks.
== configuration file settings
	-c <string>                      : Configuration file path.
	-env <string>                    : Environment file path.
== parallelization settings
	-no_par <bool>                   : Serialize all tasks (individual tasks can still use multiple threads up to '-nth').
	-nth <int>                       : Maximum # threads for a pipeline. (default: 8).
== cluster/system/resource settings
	-wt <string>                     : Walltime for all single-threaded tasks (example: 8:10:00, 3h, 3600, default: 5h50m, 5:50:00).
	-memory <string>                 : Maximum memory for all single-threaded tasks (equivalent to '-mem', example: 4.5G, 1024M, default: 7G).
	-use_system <string>             : Force to use a system (equivalent to 'bds -s [SYSTEM_NAME] ...', any system defined in bds.config can be used).
	-nice <int>                      : Set process priority for all tasks (default: 0; -20 (highest) ~ 19 (lowest) ).
	-retrial <int>                   : # of Retrial for failed tasks (default: 0).
	-q <string>                      : Submit tasks to a specified cluster queue.
	-unlimited_mem_wt <bool>         : Use unlimited max. memory and walltime.
== shell environment settings
	-mod <string>                    : Modules separated by ; (example: "bowtie/2.2.4; bwa/0.7.7; picard-tools/1.92").
	-shcmd <string>                  : Shell commands separated by ;. Shell var. must be written as ${VAR} not as $VAR (example: "export PATH=${PATH}:/usr/test; VAR=test").
	-addpath <string>                : Path separated by ; or : to be PREPENDED to \$PATH (example: "/bin/test:${HOME}/utils").
	-conda_env <string>              : Anaconda Python (or Miniconda) environment name for all softwares including Python2.
	-conda_env_py3 <string>          : Anaconda Python (or Miniconda) environment name for Python3.
	-conda_bin_dir <string>          : Anaconda Python (or Miniconda) bin directory.
	-cluster_task_min_len <int>      : Minimum length for a cluster job in seconds (dealing with NFS delayed write, default: 60).
	-cluster_task_delay <int>        : Constant delay for every job in seconds (dealing with NFS delayed write, default: 0).
== output/title settings
	-out_dir <string>                : Output directory (default: out).
	-title <string>                  : Prefix for HTML report and outputs without given prefix.
== species settings
	-species <string>                : Species. need to specify '-species_file' too if you have not installed genome database with 'install_genome_data.sh'.
	-species_file <string>           : Species file path.
	-species_browser <string>        : Species name in WashU genome browser.
	-ref_fa <string>                 : Reference genome sequence fasta.
	-chrsz <string>                  : Chromosome sizes file path (use fetchChromSizes from UCSC tools).
	-blacklist <string>              : Blacklist bed.
	-seq_dir <string>                : Reference genome sequence directory path (where chr*.fa exist).
== ENCODE accession settings
	-ENCODE_accession <string>       : ENCODE experiment accession ID (or dataset).
	-ENCODE_award_rfa <string>       : ENCODE award RFA (e.g. ENCODE3).
	-ENCODE_assay_category <string>  : ENCODE assay category.
	-ENCODE_assay_title <string>     : ENCODE assay title.
	-ENCODE_award <string>           : ENCODE award (e.g. /awards/U41HG007000/).
	-ENCODE_lab <string>             : Lab (e.g. /labs/anshul-kundaje/)
	-ENCODE_assembly <string>        : hg19, GRCh38, mm9, mm10.
	-ENCODE_alias_prefix <string>    : Alias prefix, Alias = alias_prefix: + filename + alias_suffix
	-ENCODE_alias_suffix <string>    : Alias suffix, Alias = alias_prefix: + filename + alias_suffix
== report settings
	-url_base <string>               : URL base for output directory.
	-viz_genome_coord <string>       : WashU genome browser genome coordinate (e.g. chr7:27117661-27153380).
== fastq input definition :
        Single-ended : For replicate '-fastq[REP_ID]', For control '-ctl_fastq[REP_ID]'
        Paired end : For replicate '-fastq[REP_ID]_[PAIR_ID]', For control '-ctl_fastq[REP_ID]_[PAIR_ID]'
== bam input (raw or filtered) definition :
        Raw bam : For replicate '-bam[REP_ID]', For control '-ctl_bam[REP_ID]'.
        Filtered bam : For replicate '-filt_bam[REP_ID]', For control '-ctl_filt_bam[REP_ID]'.
== tagalign input definition :
        For replicate '-tag[REP_ID]', For control '-ctl_tag[REP_ID]'.
== narrow peak input definition : 
        For true replicates, use '-peak1' and '-peak2',
        For pooled replicates, use '-peak_pooled',
        For two PR (self-pseudo-replicates), use '-peak[REP_ID]_pr1' and '-peak[REP_ID]_pr2'
        For two PPR (pooled pseudo-replicates), use '-peak_ppr1' and '-peak_ppr2'
== input endedness settings (SE or PE) :
	-se <bool>                       : Singled-ended data set. To specify it for each replicate, '-se[REP_ID]' for exp. reps, '-ctl_se[CTL_ID]' for control.
	-pe <bool>                       : Paired end data set. To specify it for each replicate, '-pe[REP_ID]' for exp. reps, '-ctl_pe[CTL_ID]' for controls.
== adapter sequence definition :
        Single-ended : For replicate '-adapter[REP_ID]'
        Paired end : For replicate '-adapter[REP_ID]_[PAIR_ID]'
== align multimapping settings
	-multimapping <int>              : # alignments reported for multimapping (default: 0).
== align bowtie2 settings (requirements: -bwt2_idx)
	-bwt2_idx <string>               : Bowtie2 index (full path prefix of *.1.bt2 file).
	-scoremin_bwt2 <string>          : Replacement --score-min for bowtie2.
	-wt_bwt2 <string>                : Walltime for bowtie2 (default: 47h, 47:00:00).
	-mem_bwt2 <string>               : Max. memory for bowtie2 (default: 12G).
	-extra_param_bwt2 <string>       : Extra parameter for bowtie2.
== adapter trimmer settings
	-adapter_err_rate <string>       : Maximum allowed adapter error rate (# errors divided by the length of the matching adapter region, default: 0.10).
	-min_trim_len <int>              : Minimum trim length for cutadapt -m, throwing away processed reads shorter than this (default: 5).
	-wt_trim <string>                : Walltime for adapter trimming (default: 23h, 23:00:00).
	-mem_trim <string>               : Max. memory for adapter trimming (default: 12G).
== postalign bam settings
	-mapq_thresh <int>               : Threshold for low MAPQ reads removal (default: 30).
	-rm_chr_from_tag <string>        : Perl style reg-ex to exclude reads from tag-aligns. (example: 'other|ribo|mito|_', '_', default: blank)
	-no_dup_removal <bool>           : No dupe removal when filtering raw bam.
	-wt_dedup <string>               : Walltime for post-alignment filtering (default: 23h, 24:00:00).
	-mem_dedup <string>              : Max. memory for post-alignment filtering (default: 12G).
	-dup_marker <string>             : Dup marker for filtering mapped reaads in BAMs: picard or sambamba (default: picard).
	-use_sambamba_markdup <bool>     : Use sambamba markdup instead of Picard MarkDuplicates (default: false).
== postalign bed/tagalign settings
	-mem_shuf <string>               : Max. memory for UNIX shuf (default: 12G).
	-fraglen0 <bool>                 : (LEGACY PARAM) Set predefined fragment length as zero for cross corr. analysis (add -speak=0 to run_spp.R).
	-speak_xcor <int>                : Set user-defined cross-corr. peak strandshift (-speak= in run_spp.R). Use -1 to disable (default: -1).
	-extra_param_xcor <string>       : Set extra parameters for run_spp.R (cross-corr. analysis only).
== postalign bed/tagalign settings
	-mem_shuf <string>               : Max. memory for UNIX shuf (default: 12G).
	-fraglen0 <bool>                 : (LEGACY PARAM) Set predefined fragment length as zero for cross corr. analysis (add -speak=0 to run_spp.R).
	-speak_xcor <int>                : Set user-defined cross-corr. peak strandshift (-speak= in run_spp.R). Use -1 to disable (default: -1).
	-extra_param_xcor <string>       : Set extra parameters for run_spp.R (cross-corr. analysis only).
== callpeak macs2 settings (requirements: -chrsz -gensz)
	-gensz <string>                  : Genome size; hs for human, mm for mouse.
	-wt_macs2 <string>               : Walltime for MACS2 (default: 23h, 23:00:00).
	-mem_macs2 <string>              : Max. memory for MACS2 (default: 15G).
	-extra_param_macs2 <string>      : Extra parameters for macs2 callpeak.
== callpeak naive overlap settings
	-nonamecheck <bool>              : bedtools intersect -nonamecheck (bedtools>=2.24.0, use this if you get bedtools intersect naming convenction warnings/errors).
== callpeak etc settings
	-npeak_filt <int>                : # top peaks filtered from a narrow peak files (default: 500000).
== IDR settings
	-idr_suffix <bool>               : Append IDR threshold to IDR output directory.
== ATAQC settings
	-tss_enrich <string>             : TSS enrichment bed for ataqc.
	-dnase <string>                  : DNase bed (open chromatin region file) for ataqc.
	-prom <string>                   : Promoter bed (promoter region file) for ataqc.
	-enh <string>                    : Enhancer bed (enhancer region file) for ataqc.
	-reg2map <string>                : Reg2map (file with cell type signals) for ataqc.
	-reg2map_bed <string>            : Reg2map_bed (file of regions used to generate reg2map signals) for ataqc.
	-roadmap_meta <string>           : Roadmap metadata for ataqc.
	-mem_ataqc <string>              : Max. memory for ATAQC (default: 20G).
	-wt_ataqc <string>               : Walltime for ATAQC (default: 47h, 47:00:00).

```
