#!/usr/bin/env nextflow

params.accessions = "assets/accessions.tsv"
params.outdir = "results"

process fondueDownload {
    publishDir "${params.outdir}/fondue", mode: 'copy'

    container 'quay.io/qiime2/pathogenome:2026.4'

    input:
        path accessions

    output:
        path "seqs.qza"
        path "metadata.tsv"

    script:
        """
        qiime fondue get-all \
            --i-accession-ids ${accessions} \
            --output-dir fondue-output
        """

}

workflow {
    accessions_ch = Channel.fromPath(params.accessions, checkIfExists: true)

    fondueDownload(accessions_ch)
}
