#!/bin/sh
# Using deepTools
if [ ! -e ../../20_piccard/Nextera-100ng-1.bam.bai ]
then
    samtools index ../../20_piccard/Nextera-100ng-1.bam
fi
alignmentSieve --maxFragmentLength 35 -b ../../20_piccard/Nextera-100ng-1.bam -o data/smallInserts.bam
samtools index data/smallInserts.bam
