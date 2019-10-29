params.sampleList = '../sampleList.txt'

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

    memory '8000 MB'

    publishDir params.BWAmemFolder, mode: 'link', overwrite: false
    
    input:
    val genome from params.genomeName
    file genome_index from bwa_index
    set sampleName, sampleId, readGroup, file(read1), file(read2) from samples_ch
    
    output:
    set sampleName, sampleId, file("${sampleId}.mappedReads.sorted.bam") into GATK_markDuplicates
    set sampleId, file("${sampleId}_BWAmem.log") into BWAmem_log
    set sampleId, file("${sampleId}_BWAmem.sh") into BWAmem_sh
   
    script:
    """
    bwa mem -t $task.cpus -R ${readGroup} ${genome} ${read1} ${read2} | \
            samtools view -h -F 4 - | \
            samtools sort -@ $task.cpus -T ./${sampleId}_tmp -O bam -o ${sampleId}.mappedReads.sorted.bam -
    cp .command.log ${sampleId}_BWAmem.log
    cp .command.sh ${sampleId}_BWAmem.sh
    """
}


process markdup {
    tag "$sampleName"

    memory '24 GB'

    publishDir params.picardFolder, mode: 'link', overwrite: false
    echo true

    input:
    set sampleName, sampleId, file(sampleId_BWA_bam) from GATK_markDuplicates.groupTuple(by:0, sort:true)

    output:
    set sampleName, file("${sampleName}.bam") into bam_for_metrics
    set sampleName, file("${sampleName}.bam") into bam_for_qualimap
    set file("${sampleName}.*Metrics*"), file("${sampleName}_markdup.log"), file("${sampleName}_markdup.sh")
 
    
    script:
    """
    ${params.javaLocation} -Djava.io.tmpdir='.' -XX:ParallelGCThreads=${task.cpus} -jar ${params.picardLocation} MarkDuplicates INPUT=${sampleId_BWA_bam.join(' INPUT=')}  OUTPUT=${sampleName}.bam METRICS_FILE=${sampleName}.MarkDuplicatesMetrics.txt TAGGING_POLICY=OpticalOnly OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 
    cp .command.log ${sampleName}_markdup.log
    cp .command.sh ${sampleName}_markdup.sh
    """
}

process metrics {
    tag "$sampleName"

    memory '8000 MB'

    publishDir params.picardFolder, mode: 'link', overwrite: false
    echo true

    input:
    set sampleName, file(bam) from bam_for_metrics


    output:
    set file("${sampleName}.*Metrics*"), file("${sampleName}_metrics.log"), file("${sampleName}_metrics.sh")
 
    
    script:
    """
    ln -s  ${params.refFolder}${params.genomeVersion}/${params.genomeName.replaceAll(/\.fn?a$/, "")}/${params.genomeName.replaceAll(/\.fn?a$/, "*")} .
    ${params.javaLocation} -Djava.io.tmpdir='.' -XX:ParallelGCThreads=${task.cpus} -jar ${params.picardLocation} CollectAlignmentSummaryMetrics REFERENCE_SEQUENCE=${params.genomeName} INPUT=$bam OUTPUT=${sampleName}.AlignmentSummaryMetrics.txt
    ${params.javaLocation} -Djava.io.tmpdir='.' -XX:ParallelGCThreads=${task.cpus} -jar ${params.picardLocation} CollectWgsMetrics REFERENCE_SEQUENCE=${params.genomeName} INPUT=$bam OUTPUT=${sampleName}.WgsMetrics.txt 
    ${params.javaLocation} -Djava.io.tmpdir='.' -XX:ParallelGCThreads=${task.cpus} -jar ${params.picardLocation} CollectInsertSizeMetrics INPUT=$bam OUTPUT=${sampleName}.InsertSizeMetrics.txt HISTOGRAM_FILE=${sampleName}.InsertSizeMetrics-Histogram.pdf 
    cp .command.log ${sampleName}_metrics.log
    cp .command.sh ${sampleName}_metrics.sh
    """
}


process qualimap {
    tag "$sampleName"

    publishDir "${params.picardFolder}/qualimap/", mode: 'link', overwrite: true
    
    input:
    set sampleName, file(bam) from bam_for_qualimap

    output:
    set sampleName, file("${sampleName}") into qm_stats
    set sampleName, file("${sampleName}_qm.log") into qm_log
    set sampleName, file("${sampleName}_qm.sh") into qm_sh
    
    script:
    """
    $params.qualimapLocation bamqc -nt $task.cpus -bam ${bam} -outdir ${sampleName} -outformat PDF:HTML --java-mem-size=15G
    cp .command.log ${sampleName}_qm.log
    cp .command.sh ${sampleName}_qm.sh
    """
}
