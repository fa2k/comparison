#!/bin/bash

# ** plotCorrelation **

# DeepTools plot correlation heatmap between samples

# Run from root analysis directory (call as plotting/plotCorrelation.sh).
# A file heatmap.pdf is created.

sample_names=`cut -f1 -d, downsample_factors.csv | tail -n +2 | sort`

echo plotCorrelation  -in 50_deeptools_summary/summary_aa_all.npz \
                --whatToPlot heatmap --corMethod spearman \
                --labels $sample_names \
                -o 50_deeptools_summary/heatmap_aa_all.pdf

# Plotcorrelation is also included in deeptools-summary.sh. Here it plots
# separately for 100ng and 10ng. The summaries are generated again in that
# script (can  be run on loki), and it outputs plotting commands similar to
# the one above, to stdout.

# ** plotCoverage **
bam_files=""
bam_names=""
for bamfile in `echo 30_downsample/*-100ng-*-aa_DS_MD.bam | sort`
do
    bam_files="$bam_files $bamfile"
    b1="${bamfile#30_downsample/}"
    b2="${b1%_aa_DS_MD.bam}"
    bam_names="$bam_names $b2"
done
echo plotCoverage -b $bam_files --label $bam_names \
                --plotFile 50_deeptools_summary/coverage_all_aa_100ng_plot.pdf \
                --outRawCounts 50_deeptools_summary/coverage_all_aa_100ng.txt

echo bamPEFragmentSize --bamfiles $bam_files \
                --histogram 50_deeptools_summary/insertSize_histogram_all_aa_100ng.pdf \
                --outRawFragmentLengths 50_deeptools_summary/insertSize_all_aa_100ng.txt

# ** plotCoverage -- 10ng **
bam_files=""
bam_names=""
for bamfile in `echo 30_downsample/*-10ng-*-aa_DS_MD.bam | sort`
do
    bam_files="$bam_files $bamfile"
    b1="${bamfile#30_downsample/}"
    b2="${b1%-aa_DS_MD.bam}"
    bam_names="$bam_names $b2"
done
echo plotCoverage -b $bam_files --label $bam_names \
                --plotFile 50_deeptools_summary/coverage_all_aa_10ng_plot.pdf \
                --outRawCounts 50_deeptools_summary/coverage_all_aa_10ng.txt

echo bamPEFragmentSize --bamfiles $bam_files \
                --histogram 50_deeptools_summary/insertSize_histogram_all_aa_10ng.pdf \
                --outRawFragmentLengths 50_deeptools_summary/insertSize_all_aa_10ng.txt

