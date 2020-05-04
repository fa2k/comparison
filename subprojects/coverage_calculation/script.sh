#!/bin/bash

REF="/ypool/bulk/fa2k/nsc/ref/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
FQ_R1="/archive/bulk/nsc/KitComparison/fastq/13-Kapa-10ng-1_S13_L005_R1_001.fastq.gz"
FQ_R2="/archive/bulk/nsc/KitComparison/fastq/13-Kapa-10ng-1_S13_L005_R2_001.fastq.gz"

PICARD="/dangerzone/applications/picard/2.22.3/picard.jar"

# Skip lines and take N_TAKE
N_SKIP=10000000
N_TAKE=100000

# Cleanup (optional)
rm -f test/*

bwa mem $REF \
        <(gunzip -c $FQ_R1 | head -n $N_SKIP | tail -n $N_TAKE) \
        <(gunzip -c $FQ_R2 | head -n $N_SKIP | tail -n $N_TAKE) \
        > test/sample.sam

samtools sort -o test/sample.sorted.bam test/sample.sam
samtools index test/sample.sorted.bam

java -jar $PICARD CollectWgsMetrics \
        R=$REF \
        I=test/sample.sorted.bam \
        O=test/sample_wgs_metrics.txt 

java -jar $PICARD CollectInsertSizeMetrics \
        I=test/sample.sorted.bam \
        H=test/sample_insert_size_hist.pdf \
        O=test/sample_insert_size_metrics.txt

java -jar $PICARD MarkDuplicates \
        I=test/sample.sorted.bam \
        O=test/sample.markdup.bam \
        METRICS_FILE=test/MarkDuplicatesMetrics.txt \
        TAGGING_POLICY=OpticalOnly \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500

samtools index test/sample.markdup.bam

java -jar $PICARD CollectWgsMetrics \
        R=$REF \
        I=test/sample.markdup.bam \
        O=test/sample_wgs_metrics_post_markdup.txt 


# Remove the first 14 alignments, but keep the header
samtools view -H test/sample.markdup.bam > test/tail.sam
samtools view test/sample.markdup.bam | tail -n+15 >> test/tail.sam

java -jar $PICARD CollectWgsMetrics \
        R=$REF \
        I=test/tail.sam \
        O=test/tail_wgs_metrics.txt
