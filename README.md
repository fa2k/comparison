## Analysis overview ##


# 1. main.nf: Alignment and pre-processing

The input is raw FASTQ data, and this is handled in main.nf. All files from each sample are
processed together in the MarkDuplicates run, so they are combined into one file per sample.


# 2. main_b.nf: Downsampling and initial analysis

Before running the next workflow, the downsampling factors must be computed. For deduplicated
data, the script compute_downsample_deduped.py computes the values and writes directly to
the expected location downsample_factors_dd_fix.txt. For non-deduplicated reads, the
normalisation is based on the number of reads in the fastq fiels. The downsample factors
file is generated using the notebook in 01_data_qc (TODO: remove "deduplicated" branch, not
used)

The number of reads is reduced to a fixed value for all samples (downsampling), using a
Picard tool. MarkDuplicates is repeated on this data, and then variant calling is performed
using GATK/HaplotypeCaller. In parallel, text files are created showing the number of reads
in the BAM files, in categories of unmapped, multi-mapped and uniquely mapped.


# 3. main_c_variants.nf: Variant call analysis

Runs a tool to score and filter variants: CNNScoreVariants, and then hap.py to determine the
concordance with the well known truth for this sample.


# 4. main_d_deeptools.nf: Deeptools

Creates Bigwig files based on the BAMs, and computes other stats from deepTools.

