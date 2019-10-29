
bam_files=""
bam_names=""
for bamfile in `echo 40_bigwig/*-100ng-*-aa.bw | sort`
do
    bam_files="$bam_files $bamfile"
    b1="${bamfile#40_bigwig/}"
    b2="${b1%.bw}"
    bam_names="$bam_names $b2"
done

scl enable rh-python35 "multiBigwigSummary bins -p 16 -b $bam_files -o 50_deeptools_summary/summary_aa_all_100ng.npz"

echo plotCorrelation  -in 50_deeptools_summary/summary_aa_all_100ng.npz \
                --whatToPlot heatmap --corMethod spearman \
                --labels $bam_names -o 50_deeptools_summary/heatmap_aa_all_100ng.pdf

bam_files=""
bam_names=""
for bamfile in `echo  40_bigwig/*-10ng-*-aa.bw | sort`
do
    bam_files="$bam_files $bamfile"
    b1="${bamfile#40_bigwig/}"
    b2="${b1%.bw}"
    bam_names="$bam_names $b2"
done

scl enable rh-python35 "multiBigwigSummary bins -p 16 -b $bam_files -o 50_deeptools_summary/summary_aa_all_10ng.npz"

echo plotCorrelation  -in 50_deeptools_summary/summary_aa_all_10ng.npz \
                --whatToPlot heatmap --corMethod spearman \
                --labels $bam_names -o 50_deeptools_summary/heatmap_aa_all_10ng.pdf
