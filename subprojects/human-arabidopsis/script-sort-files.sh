#!/bin/sh

mkdir -p data tmp

SAMPLE=NEB-100ng-1

singularity run -B ../../20_piccard/$SAMPLE.bam:/input.bam \
                -B $PWD/data/:/data \
                ../../container-images/picard-2.20.4.simg SortSam \
                I=/input.bam \
                O=/data/$SAMPLE-human-qsort.bam \
                SO=queryname \
                TMP_DIR=tmp/ \
                CREATE_INDEX=true


singularity run -B ../../arabidopsis/20_piccard/$SAMPLE.bam:/input.bam \
                -B $PWD/data/:/data \
                ../../container-images/picard-2.20.4.simg SortSam \
                I=/input.bam \
                O=/data/$SAMPLE-arabidopsis-qsort.bam \
                SO=queryname \
                TMP_DIR=tmp/ \
                CREATE_INDEX=true
