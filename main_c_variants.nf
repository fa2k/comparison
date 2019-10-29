
Channel
    .fromPath('downsample_factors_fix.txt')
    .splitCsv(header:true, sep:"\t")
    .map{ row -> tuple(
                    row.sampleName,
                    file("${params.vcfFolder}/${row.sampleName}.vcf"),
                    file("${params.downSampleFolder}/${row.sampleName}_DS_MD.bam"),
                    file("${params.downSampleFolder}/${row.sampleName}_DS_MD.bam.bai")
                ) }
    .set { samples_ch }


process score {
    tag "$sampleName"

    clusterOptions "-C el7,avx"
    cpus 8
    memory '15 GB'
    time '12d'

    publishDir "${params.vcfFolder}", mode: 'link', overwrite: true

    input:
    set sampleName, file(vcf), file(bam), file(bai) from samples_ch

    output:
    set sampleName, file("${sampleName}.scored.vcf") into scored_vcfs
    file("${sampleName}.scored.vcf.idx")

    script:
    """
    ln -s ${params.refFolder}${params.genomeVersion}/${params.genomeName.replaceAll(/.fn?a$/, "")}/${params.genomeName.replaceAll(/\.fn?a$/, "*")} .
    # It's impossible to resume this one in nextflow, so we skip the ones that are already done
    if [[ -e "${params.vcfFolder}/${sampleName}.scored.vcf" ]]
    then
        cp "${params.vcfFolder}/${sampleName}.scored.vcf" .
        cp "${params.vcfFolder}/${sampleName}.scored.vcf.idx" .
    else
        gatk CNNScoreVariants \
           --input $bam \
           --variant $vcf \
           --reference $params.genomeName \
           --output ${sampleName}.scored.vcf \
           --tensor-type read_tensor
    fi

    cp .command.log ${sampleName}_rc.log
    cp .command.sh ${sampleName}_rc.sh
    """
}

process happy {
    tag "$sampleName"

    clusterOptions "-C el7"
    beforeScript """ln ${params.refFolder}${params.genomeVersion}/${params.genomeName.replaceAll(/.fn?a$/, "")}/${params.genomeName.replaceAll(/\.fn?a$/, "*")} .; ln ${params.happyConfidentCallsBasePath}* ."""

    publishDir "${params.variantAnalysisFolder}", mode: 'link', overwrite: true

    input:
    set sampleName, file(vcf) from scored_vcfs

    output:
    file("${sampleName}_happy")


    script:
    """
    mkdir tmp ${sampleName}_happy
    export TMPDIR=./tmp
    /opt/hap.py/bin/hap.py --threads $task.cpus -r $params.genomeName -f $params.happyConfidentCallsBed \
        --roc CNN_2D $params.happyConfidentCallsVcf $vcf -o ${sampleName}_happy/${sampleName}

    cp .command.log ${sampleName}_hap.log
    cp .command.sh ${sampleName}_hap.sh
    """
}

