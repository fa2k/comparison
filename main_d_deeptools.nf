Channel
    .fromPath('downsample_factors_fix.txt')
    .splitCsv(header:true, sep:"\t")
    .map{ row -> tuple(
                    "aa",
                    row.sampleName,
                    file("${params.downSampleFolder}/${row.sampleName}_DS_MD.bam"),
                    file("${params.downSampleFolder}/${row.sampleName}_DS_MD.bam.bai")
                ) }
    .into { samples1; samples_ch_gc }

Channel
    .fromPath('downsample_factors_dd_fix.txt')
    .splitCsv(header:true, sep:"\t")
    .map{ row -> tuple(
                    "dd",
                    row.sampleName,
                    file("${params.downSampleFolder}/${row.sampleName}_DS_MD.bam"),
                    file("${params.downSampleFolder}/${row.sampleName}_DS_MD.bam.bai")
                ) }
    .set { samples2 }

samples1.concat(samples2).into { samples_ch_bw1; samples_ch_bw2 }


//process bigwig {
//    tag "$sampleName"
//
//    cpus 8
//    memory '16 GB'
//    clusterOptions '-C rhel'
//
//    publishDir "${params.bigwigFolder}", mode: 'link', overwrite: true
//
//    input:
//    set dedup, sampleName, file(bam), file(bai) from samples_ch_bw1
//
//    output:
//    set dedup, val("all"), file("${sampleName}.bw") into multibwsummary_all
//
//    script:
//    """
//    scl enable rh-python35 "bamCoverage \
//                        --binSize 1 -p $task.cpus \
//                        -b $bam -o ${sampleName}.bw"
//
//    cp .command.log ${sampleName}_bw.log
//    cp .command.sh ${sampleName}_bw.sh
//    """
//}
//
//
//process bigwigq30 {
//    tag "$sampleName"
//
//    cpus 16
//    memory '16 GB'
//    clusterOptions "-C rhel"
//
//    publishDir "${params.bigwigFolder}", mode: 'link', overwrite: true
//
//    input:
//    set dedup, sampleName, file(bam), file(bai) from samples_ch_bw2
//
//    output:
//    set dedup, val("q30"), file("${sampleName}_mapQ30.bw") into multibwsummary_q30
//
//    script:
//    """
//    scl enable rh-python35 "bamCoverage \
//                        --binSize 1 -p $task.cpus --minMappingQuality 30 \
//                        -b $bam -o ${sampleName}_mapQ30.bw"
//
//    cp .command.log ${sampleName}_bwq30.log
//    cp .command.sh ${sampleName}_bwq30.sh
//    """
//}
//
//multibwsummaries = multibwsummary_all.concat(multibwsummary_q30)
//
//process bwsummary {
//    cpus 16
//    memory '16 GB'
//    clusterOptions "-C rhel"
//
//    publishDir params.deeptoolsSummaryFolder, mode: 'link', overwrite: true
//
//    input:
//    set dedup, type, file(bigwig) from multibwsummaries.groupTuple(by: [0,1], sort: true)
//
//    output:
//    file("summary_${dedup}_${type}.npz") into bwsummary_out
//
//    script:
//    """
//    scl enable rh-python35 "multiBigwigSummary bins \
//                        -p $task.cpus -b $bigwig \
//                        -o summary_${dedup}_${type}.npz"
//
//    cp .command.log bws_${dedup}_${type}.log
//    cp .command.sh bws_${dedup}_${type}.sh
//    """
//}

process gcbias {
    tag "$sampleName"

    cpus 4
    memory '15 GB'
    clusterOptions "-C rhel"

    publishDir "${params.gcBiasFolder}", mode: 'link', overwrite: true

    input:
    set dedup, sampleName, file(bam), file(bai) from samples_ch_gc

    output:
    file("${sampleName}.txt")

    script:
    """
    /usr/bin/scl enable rh-python36 "computeGCBias -g $params.refFolder/2bit/hg38.2bit \
                            --effectiveGenomeSize 2913022398 \
                            -p $task.cpus -b $bam -o ${sampleName}.txt"

    cp .command.log ${sampleName}_gc.log
    cp .command.sh ${sampleName}_gc.sh
    """
}

// plotCorrelation handled in separate script, see plotting/deeptools.sh
// It doesn't work on SCL python35, though it would have worked on python36,
// added to support computeGCBias.


//2913022398
