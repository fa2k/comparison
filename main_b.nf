
samples1 = Channel
    .fromPath('downsample_factors.txt')
    .splitCsv(header:true, sep:"\t")
    .map{ row -> tuple(row.sampleName, file(row.fileBAM), row.downSampleFactor) }

//// Deduplicated data -- this section is disabled because it is not used at the moment.
//// Instead: process samples1 as samples_ch

samples1.set { samples_ch }

// samples2 = Channel
//     .fromPath('downsample_factors_dd_fix.txt')
//     .splitCsv(header:true, sep:"\t")
//     .map{ row -> tuple(row.sampleName, file(row.fileBAM), row.downSampleFactor) }
// 
// samples1.concat(samples2).set { samples_ch }


process downsample {
    tag "$sampleName"

    // Does not publish the .bam, which is not needed (we publish the MarkDuplicates bam instead)
    publishDir params.downSampleFolder, mode: 'link', overwrite: true, pattern: '*.txt'

    input:
    set sampleName, file(bam), downSampleFactor from samples_ch
    
    output:
    set sampleName, file("${sampleName}_DS.bam") into bam_for_markdup
    set file("${sampleName}_*"), file("${sampleName}.DownSampleSam.txt")
    
    script:
    """
    ${params.javaLocation} -Djava.io.tmpdir='.' -XX:ParallelGCThreads=${task.cpus} -jar ${params.picardLocation} DownsampleSam INPUT=${bam} OUTPUT=${sampleName}_DS.bam METRICS_FILE=${sampleName}.DownSampleSam.txt PROBABILITY=${downSampleFactor}

    ln -s  ${params.refFolder}${params.genomeVersion}/${params.genomeName.replaceAll(/.fn?a$/, "")}/${params.genomeName.replaceAll(/\.fn?a$/, "*")} . 
    
    cp .command.log ${sampleName}_DownSampleSam.log
    cp .command.sh ${sampleName}_DownSampleSam.sh    
    """
}

process markdup {
    tag "$sampleName"

    memory '15 GB'
    publishDir params.downSampleFolder, mode: 'link', overwrite: true

    input:
    set sampleName, file(bam) from bam_for_markdup

    output:
    set sampleName, file("${sampleName}_DS_MD.bam") into  bam_for_qualimap, bam_for_vc, bam_for_samtools, bam_for_picard, bam_for_indexing, bam_for_cnn_variantscore
    file "${sampleName}_DS_MD.bam" into count_reads_md
    file "${sampleName}_*"

    script:
    """
    ${params.javaLocation} -Djava.io.tmpdir='.' -XX:ParallelGCThreads=${task.cpus} -jar ${params.picardLocation} MarkDuplicates INPUT=$bam  OUTPUT=${sampleName}_DS_MD.bam METRICS_FILE=${sampleName}_DS.MarkDuplicatesMetrics.txt TAGGING_POLICY=OpticalOnly OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500

    cp .command.log ${sampleName}_MarkDuplicates.log
    cp .command.sh ${sampleName}_MarkDuplicates.sh    
    """
}

process samtools {
    tag "$sampleName"

    cpus 1
    memory '1 GB'

    publishDir params.downSampleFolder, mode: 'link', overwrite: true

    input:
    set sampleName, file(bam) from bam_for_samtools

    output:
    file "${sampleName}_DS_MD_mapped.bam" into count_reads_st1
    file "${sampleName}_DS_MD_mapped_multiple.bam" into count_reads_st2
    file "${sampleName}_DS_MD_mapped_single.bam" into count_reads_st3
    file "${sampleName}_samtools.*"


    script:
    """
    samtools view -b -h -F 4 $bam >  ${sampleName}_DS_MD_mapped.bam
    samtools view -b -h -q 1 ${sampleName}_DS_MD_mapped.bam -U ${sampleName}_DS_MD_mapped_multiple.bam >  ${sampleName}_DS_MD_mapped_single.bam 

    cp .command.log ${sampleName}_samtools.log
    cp .command.sh ${sampleName}_samtools.sh    
    """
}

count_reads = count_reads_md.mix(count_reads_st1, count_reads_st2, count_reads_st3)

process count_reads {
    cpus 1
    memory '1 GB'

    publishDir params.downSampleFolder, mode: 'link', overwrite: true

    input:
    file bam from count_reads

    output:
    file "${bam}.ReadCount.txt"

    script:
    """
    samtools view -c -F 0x900 $bam > ${bam}.ReadCount.txt
    """
}

process metrics {
    tag "$sampleName"

    cpus 4
    
    publishDir params.downSampleFolder, mode: 'link', overwrite: true

    input:
    set sampleName, file(bam) from bam_for_picard
    
    output:
    file "${sampleName}_*.txt"
    file "${sampleName}_*.pdf"
    file "${sampleName}_PicardMetrics.*"
    
    script:
    """
    ln -s  ${params.refFolder}${params.genomeVersion}/${params.genomeName.replaceAll(/.fn?a$/, "")}/${params.genomeName.replaceAll(/\.fn?a$/, "*")} . 
    
    ${params.javaLocation} -Djava.io.tmpdir='.' -XX:ParallelGCThreads=${task.cpus} -jar ${params.picardLocation} CollectAlignmentSummaryMetrics REFERENCE_SEQUENCE=${params.genomeName} INPUT=$bam OUTPUT=${sampleName}_DS_MD.AlignmentSummaryMetrics.txt
    ${params.javaLocation} -Djava.io.tmpdir='.' -XX:ParallelGCThreads=${task.cpus} -jar ${params.picardLocation} CollectWgsMetrics REFERENCE_SEQUENCE=${params.genomeName} INPUT=$bam OUTPUT=${sampleName}_DS_MD.WgsMetrics.txt 
    ${params.javaLocation} -Djava.io.tmpdir='.' -XX:ParallelGCThreads=${task.cpus} -jar ${params.picardLocation} CollectInsertSizeMetrics INPUT=$bam OUTPUT=${sampleName}_DS_MD.InsertSizeMetrics.txt HISTOGRAM_FILE=${sampleName}_DS_MD.InsertSizeMetrics-Histogram.pdf 

    cp .command.log ${sampleName}_PicardMetrics.log
    cp .command.sh ${sampleName}_PicardMetrics.sh    
    """
}

process qualimap {
    tag "$sampleName"

    publishDir "${params.downSampleFolder}/qualimap/", mode: 'link', overwrite: true
    
    input:
    set sampleName, file(bam) from bam_for_qualimap

    output:
    set sampleName, file("${sampleName}_DS") into qm_stats
    set sampleName, file("${sampleName}_qm.log") into qm_log
    set sampleName, file("${sampleName}_qm.sh") into qm_sh
    
    script:
    """
    $params.qualimapLocation bamqc -nt $task.cpus -bam ${bam} -outdir ${sampleName}_DS -outformat PDF:HTML --java-mem-size=15G
    cp .command.log ${sampleName}_qm.log
    cp .command.sh ${sampleName}_qm.sh
    """
}

process indexing {
    publishDir params.downSampleFolder, mode: 'link', overwrite: true
    
    cpus 4
    memory '1 GB'

    input:
    set sampleName, file(bam) from bam_for_indexing

    output:
    set sampleName, file("${bam}.bai") into bam_index
    
    script:
    """
    samtools index $bam
    """
}

process hc {
    tag "$sampleName"

    time '7d'

    publishDir "${params.vcfFolder}", mode: 'link', overwrite: true
    
    input:
    set sampleName, file(bam), file(bai) from bam_for_vc.join(bam_index) // Join on sample name

    output:
    set sampleName, file("${sampleName}.vcf") into hc_vcfs
    file "${sampleName}_hc.log"
    file "${sampleName}_hc.sh"
    
    script:
    """
    ln -s  ${params.refFolder}${params.genomeVersion}/${params.genomeName.replaceAll(/.fn?a$/, "")}/${params.genomeName.replaceAll(/\.fn?a$/, "*")} . 
    export JAVA_HOME=/data/common/tools/jdk8/current
    export PATH=/data/common/tools/jdk8/current/bin:$PATH
    $params.gatkLocation --java-options "-Xmx4g" HaplotypeCaller \
                -R $params.genomeName \
                -I $bam \
                -O ${sampleName}.vcf
    cp .command.log ${sampleName}_hc.log
    cp .command.sh ${sampleName}_hc.sh
    """
}

//process score {
//    tag "$sampleName"
//
//
//    publishDir "${params.vcfFolder}", mode: 'link', overwrite: true
//
//    input:
//    set sampleName, file(vcf) from hc_vcfs
//    file bam_for_cnn_variantscore
//
//    output:
//    set sampleName, file("${sampleName}.scored.vcf}") into scored_vcfs
//	
//
//    script:
//    """
//    ln -s ${params.refFolder}${params.genomeVersion}/${params.genomeName.replaceAll(/.fn?a$/, "")}/${params.genomeName.replaceAll(/\.fn?a$/, "*")} . 
//    $params.gatkLocation CNNScoreVariants \
//       -I $bam_for_cnn_variantscore \
//       -V $vcf \
//       -R $params.genomeName \
//       -O ${sampleName}.scored.vcf \
//       -tensor-type read-tensor
//
//       #-inference-batch-size 2 \
//       #-transfer-batch-size 2 \
//
//    cp .command.log ${sampleName}_rc.log
//    cp .command.sh ${sampleName}_rc.sh
//    """
//}


