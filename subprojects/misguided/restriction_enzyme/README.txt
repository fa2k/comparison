Determine if Quanta kit has higher duplicates due to a non-random restriction enzyme.

Plan:
    - Compare nucleotide composition at start of R1 and R2 in non-optical duplicates, and in all reads.
    - Compare with other kit.


# Non-sequencing duplicates were extraced using this command

samtools view -h Quanta-100ng-4.bam -f 0x400 | grep -vw 'DT:Z:SQ' | samtools view -b - > Quanta-100ng-4.duplicates.nonSQ.bam


# Kit for comparison, with lower duplicates: Swift2S.
# Swift has a clear difference in duplication with 10 and 100 ng as expected.
# The percentages are (100ng): Quanta 2.2 %, Swift ca. 1 %

samtools view -h Swift2S-100ng-1.bam -f 0x400 | grep -vw 'DT:Z:SQ' | samtools view -b - > Swift2S-100ng-1.duplicates.nonSQ.bam

# Run fastqc on all (it was done individually, that's equalivant)
fastqc *.bam

# MultiQC
multiqc bamqc -bam F.bam



