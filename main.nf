params.sampleList = 'sampleList.txt'

Channel
    .fromPath(params.sampleList)
    .splitCsv(header:true, sep:"\t")
    .map{ row -> tuple(row.sampleName, row.sampleId, row.readGroup, file(row.read1), file(row.read2)) }
    .set { samples_ch }


process bwa_buildIndex {
    tag "copying BWA MEM index files"
    
    input:
    val genome from params.genomeName
    
    output:
    file "${genome}.*" into bwa_index
    
    script:    
    """
    ln -s  ${params.refFolder}/bwa_index/${genome}* . 
    """
}


process bwa_mem {
    tag "$sampleId"
    publishDir params.BWAmemFolder, mode: 'copy', overwrite: false
    
    input:
    val genome from params.genomeName
    file genome_index from bwa_index
    set sampleName, sampleId, readGroup, file(read1), file(read2) from samples_ch
    
    output:
    set sampleName, sampleId, file("${sampleId}.sorted.bam") into GATK_markDuplicates
    set sampleId, file("${sampleId}_BWAmem.log") into BWAmem_log
    set sampleId, file("${sampleId}_BWAmem.sh") into BWAmem_sh
   
    script:
    """
    bwa mem -t $params.cpu -R ${readGroup} ${genome} ${read1} ${read2} | samtools sort -m 4G -@ $params.cpu -T ./${sampleId}_tmp -O bam -o ${sampleId}.sorted.bam -
    cp .command.log ${sampleId}_BWAmem.log
    cp .command.sh ${sampleId}_BWAmem.sh
    """
}


process piccard {
    tag "$sampleName"
    publishDir params.picardFolder, mode: 'copy', overwrite: false
    echo true

    input:
    set sampleName, sampleId, file(sampleId_BWA_bam) from GATK_markDuplicates.groupTuple(by:0, sort:true)

    output:
    set sampleName, file("${sampleName}.bam") into bam_for_qualimap
    set sampleName, file("${sampleName}.bam") into bam_for_deduplicate
    set file("${sampleName}.*Metrics*"), file("${sampleName}_picard.log"), file("${sampleName}_picard.sh")
 
    
    script:
    """
    ln -s  ${params.refFolder}${params.genomeVersion}/${params.genomeName.replaceAll(/\.fn?a$/, "")}/${params.genomeName.replaceAll(/\.fn?a$/, "*")} .
    ${params.javaLocation} -Djava.io.tmpdir='.' -XX:ParallelGCThreads=${params.cpu} -jar ${params.picardLocation} MarkDuplicates INPUT=${sampleId_BWA_bam.join(' INPUT=')}  OUTPUT=${sampleName}.bam METRICS_FILE=${sampleName}.MarkDuplicatesMetrics.txt TAGGING_POLICY=OpticalOnly OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 
    ${params.javaLocation} -Djava.io.tmpdir='.' -XX:ParallelGCThreads=${params.cpu} -jar ${params.picardLocation} CollectAlignmentSummaryMetrics REFERENCE_SEQUENCE=${params.genomeName} INPUT=${sampleName}.bam OUTPUT=${sampleName}.AlignmentSummaryMetrics.txt
    ${params.javaLocation} -Djava.io.tmpdir='.' -XX:ParallelGCThreads=${params.cpu} -jar ${params.picardLocation} CollectWgsMetrics REFERENCE_SEQUENCE=${params.genomeName} INPUT=${sampleName}.bam OUTPUT=${sampleName}.WgsMetrics.txt 
    ${params.javaLocation} -Djava.io.tmpdir='.' -XX:ParallelGCThreads=${params.cpu} -jar ${params.picardLocation} CollectInsertSizeMetrics INPUT=${sampleName}.bam OUTPUT=${sampleName}.InsertSizeMetrics.txt HISTOGRAM_FILE=${sampleName}.InsertSizeMetrics-Histogram.pdf 
    cp .command.log ${sampleName}_picard.log
    cp .command.sh ${sampleName}_picard.sh
    """
}

// exec:
// println "-I ${sampleId_BWA_bam.join(' -I ')}"


process dedup {
    tag "$sampleName"
    publishDir params.dedupFolder, mode: 'copy', overwrite: false

    echo true

    input:
    set sampleName, file(bam) from bam_for_deduplicate

    output:
    set sampleName, file("${sampleName}.filter.bam") into dedup_bam_for_counting
    set file("${sampleName}.filter.bam"), file("${sampleName}_dedup.log"), file("${sampleName}_dedup.sh")

    script:
    """
    samtools view -@${task.cpus} -h $bam | grep -vw 'DT:Z:SQ' | samtools view -@${task.cpus} -b - > ${sampleName}.filter.bam
    cp .command.log ${sampleName}_dedup.log
    cp .command.sh ${sampleName}_dedup.sh
    """
}

process readCount {
    tag "$sampleName"
    publishDir params.dedupFolder, mode: 'copy', overwrite: false

    echo true

    input:
    set sampleName, file(bam) from dedup_bam_for_counting

    output:
    file "${sampleName}.readCount.txt"

    script:
    """
    samtools view -c -F 0x900 $bam > ${sampleName}.readCount.txt
    cp .command.log ${sampleName}_dedup.log
    cp .command.sh ${sampleName}_dedup.sh
    """
}


