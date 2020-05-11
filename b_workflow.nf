//// NOTE: Uses nextflow second version language.
nextflow.preview.dsl=2

//// Workflow "b" -- downsample based on coverage.
//// Contrary to the other workflow, this one performs all of the analysis in a single
//// Nextflow file.

// Usage information:
// Inputs are BAM files from 20_piccard.
// The file b_downsample_factors.csv should contain a list of samples, paths and downsampling factors.
// The CSV file can be generated using the notebook b30_coverage_downsampling_factors.ipynb


process downsample {
    tag "$sampleName"
    label "picard_container"

    input:
    tuple sampleName, file(bam), downSampleFactor
    
    output:
    tuple sampleName, file("${sampleName}_DS.bam")
    
    script:
    """
    $params.picardCommand DownsampleSam \
            INPUT=${bam} OUTPUT=${sampleName}_DS.bam \
            METRICS_FILE=${sampleName}.DownSampleSam.txt \
            PROBABILITY=${downSampleFactor}
    """
}

process markdup {
    tag "$sampleName"
    label "picard_container"
    memory '15 GB'
    publishDir "${params.workflowBranchId}30_downsample", mode: 'link', overwrite: true

    input:
    tuple sampleName, file(bam)

    output:
    tuple sampleName, file("${sampleName}_DS_MD.bam")

    script:
    """
    $params.picardCommand MarkDuplicates \
            INPUT=$bam OUTPUT=${sampleName}_DS_MD.bam \
            METRICS_FILE=${sampleName}_DS.MarkDuplicatesMetrics.txt \
            TAGGING_POLICY=OpticalOnly \
            OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500
    """
}

process metrics {
    tag "$sampleName"
    label "picard_container"
    cpus 4
    publishDir "${params.workflowBranchId}30_downsample", mode: 'link', overwrite: true

    input:
    path 'genome.fa'
    tuple sampleName, file(bam)
    
    output:
    file "${sampleName}_*.txt"
    file "${sampleName}_*.pdf"
    
    script:
    """
    $params.picardCommand CollectAlignmentSummaryMetrics \
                REFERENCE_SEQUENCE=genome.fa \
                INPUT=$bam OUTPUT=${sampleName}_DS_MD.AlignmentSummaryMetrics.txt
    $params.picardCommand CollectWgsMetrics \
                REFERENCE_SEQUENCE=genome.fa \
                INPUT=$bam OUTPUT=${sampleName}_DS_MD.WgsMetrics.txt 
    $params.picardCommand CollectInsertSizeMetrics \
                INPUT=$bam OUTPUT=${sampleName}_DS_MD.InsertSizeMetrics.txt \
                HISTOGRAM_FILE=${sampleName}_DS_MD.InsertSizeMetrics-Histogram.pdf
    """
}

process indexing {
    publishDir "${params.workflowBranchId}30_downsample", mode: 'link', overwrite: true
    label 'samtools_container'
    cpus 4
    memory '1 GB'

    input:
    tuple sampleName, file(bam)

    output:
    tuple sampleName, file("${bam}.bai")
    
    script:
    """
    samtools index $bam
    """
}

process hc {
    tag "$sampleName"
    label 'gatk_container'
    time '7d'

    input:
    path 'genome.fa'
    path 'genome.fa.fai'
    path 'genome.dict'
    tuple sampleName, file(bam), file(bai)

    output:
    tuple sampleName, file("${sampleName}.vcf")
    
    script:
    """
    $params.gatkCommand HaplotypeCaller \
                -R genome.fa \
                -I $bam \
                -O ${sampleName}.vcf
    """
}

process score {
    tag "$sampleName"
    label 'gatk_container'
    memory '15 GB'
    time '12d'

    publishDir "${params.workflowBranchId}40_vcf", mode: 'link', overwrite: true

    input:
    path 'genome.fa'
    path 'genome.fa.fai'
    path 'genome.dict'
    tuple sampleName, file(vcf), file(bam), file(bai)

    output:
    tuple sampleName, file("${sampleName}.scored.vcf"), file("${sampleName}.scored.vcf.idx")

    script:
    """
    $params.gatkCommand CNNScoreVariants \
       -I $bam \
       -V $vcf \
       -R genome.fa \
       -O ${sampleName}.scored.vcf \
       -tensor-type read_tensor
    """
}

process happy {
    tag "$sampleName"
    label 'happy_container'

    publishDir "${params.workflowBranchId}50_variant_analysis", mode: 'link', overwrite: true

    input:
    path 'genome.fa'
    path 'genome.fa.fai'
    path 'cc.vcf'
    path 'cc.bed'
    tuple sampleName, file(vcf), file(vcfi)

    output:
    file("${sampleName}_happy")

    script:
    """
    mkdir ${sampleName}_happy
    /opt/hap.py/bin/hap.py --threads $task.cpus -r genome.fa -f cc.bed \
        --roc CNN_2D cc.vcf $vcf -o ${sampleName}_happy/${sampleName}
    """
}

process bigwigqx {
    tag "$sampleName"
    label 'deeptools_container'
    cpus 8
    memory '16 GB'

    publishDir "${params.workflowBranchId}40_bigwig", mode: 'link', overwrite: true

    input:
    tuple val(qthreshold), sampleName, file(bam), file(bai)

    output:
    tuple qthreshold, sampleName, file("${sampleName}_mapQ${qthreshold}.bw")

    script:
    """
    bamCoverage --binSize 1 -p $task.cpus --minMappingQuality $qthreshold \
                        -b $bam -o ${sampleName}_mapQ${qthreshold}.bw
    """
}

// TODO -- maybe summary file is not needed
/* process bwsummary {
    label 'deeptools_container'
    cpus 8
    memory '16 GB'

    publishDir "${params.workflowBranchId}50_deeptools_summary", mode: 'link', overwrite: true

    input:
    tuple sampleName, qthreshold, val(sampleNames), file(bigwigs)

    output:
    file("summary_${qthreshold}.npz") into bwsummary_out

    script:
    """
    multiBigwigSummary bins -p $task.cpus -b $bigwigs -o summary_${qthreshold}.npz
    """
} */

process gcbias {
    tag "$sampleName"
    label 'deeptools_container'
    cpus 4
    memory '15 GB'

    publishDir "${params.workflowBranchId}40_gc_bias", mode: 'link', overwrite: true

    input:
    path 'genome.2bit'
    tuple sampleName, file(bam), file(bai)

    output:
    file("${sampleName}.txt")

    script:
    """
    computeGCBias -g genome.2bit --effectiveGenomeSize 2913022398 \
                            -p $task.cpus -b $bam -o ${sampleName}.txt
    """
}

process plotCoverage {
    label 'deeptools_container'
    cpus 4
    memory '15 GB'

    publishDir "${params.workflowBranchId}50_deeptools_summary", mode: 'link', overwrite: true

    input:
    tuple val(sampleNames), val(conc), file(bam), file(bai)

    output:
    file("coverage_${conc}*")

    script:
    """
    plotCoverage -b $bam --label ${sampleNames.join(' ')} \
            --plotFile coverage_${conc}_plot.pdf \
            --outRawCounts coverage_${conc}.txt
    """
}

process bamPEFragmentSize {
    label 'deeptools_container'
    cpus 4
    memory '15 GB'

    publishDir "${params.workflowBranchId}50_deeptools_summary", mode: 'link', overwrite: true

    input:
    tuple val(sampleNames), val(conc), file(bam), file(bai)

    output:
    file("insertSize_${conc}*")

    script:
    """
    bamPEFragmentSize --bamfiles $bam -samplesLabel ${sampleNames.join(' ')} \
                --numberOfProcessors 4 --binSize 3000 --distanceBetweenBins 10000 \
                --histogram insertSize_${conc}_plot.pdf \
                --outRawFragmentLengths insertSize_${conc}.txt
    """
}

genome = Channel.value(params.refPath)
genomeFai = Channel.value("${params.refPath}.fai")
genomeDict = Channel.value(params.refPath.replaceAll(/\.fn?a$/, '.dict'))
genome2bit = Channel.value(params.ref2bitPath)
confidentCallsBed = Channel.value(params.confidentCallsBedPath)
confidentCallsVcf = Channel.value(params.confidentCallsVcfPath)
confidentCallsIndexes = confidentCallsBed.concat(confidentCallsVcf).collect { "${it}.tbi"}

qThresholds = Channel.fromList([0, 20])

inputSamples = Channel
    .fromPath('b_downsample_factors.csv')
    .splitCsv(header:true)
    .map{ row -> tuple(row.LIBRARY, file(row.Path), row.Scaling) }

concentrations = inputSamples.map { it[0] }.splitCsv(sep: "-").map { it[1] }
libraryConc = inputSamples.map { it[0] }.merge(concentrations)

workflow {
    downsample(inputSamples)
    markdup(downsample.out)
    metrics(genome, markdup.out)
    indexing(markdup.out)
    bamWithIndex = markdup.out.join(indexing.out)
    hc(genome, genomeFai, genomeDict, bamWithIndex)
    score(genome, genomeFai, genomeDict, hc.out.join(bamWithIndex))
    happy(genome, genomeFai, confidentCallsVcf, confidentCallsBed, score.out)
    bigwigqx(qThresholds.combine(bamWithIndex))
    //TODO or not TODO: multibwsummaries.groupTuple(by: [0,1], sort: true)
    bamWithConcAndIndex = libraryConc.join(bamWithIndex)
    plotCoverage(bamWithConcAndIndex.groupTuple(by: 1))
    gcbias(genome2bit, bamWithIndex)
}