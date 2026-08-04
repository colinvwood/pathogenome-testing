#!/usr/bin/env nextflow

params.accessions = null
params.fondue_email = null

params.import_threads = 1
params.import_memory_gb = 8
params.import_time_hours = 1
params.fondue_threads = 8
params.fondue_memory_gb = 16
params.fondue_time_hours = 36
params.megahit_threads = 12
params.megahit_memory_gb = 64
params.megahit_time_hours = 48
params.prodigal_threads = 12
params.prodigal_memory_gb = 32
params.prodigal_time_hours = 24
params.amr_db_threads = 8
params.amr_db_memory_gb = 16
params.amr_db_time_hours = 24
params.amr_threads = 4
params.amr_memory_gb = 128
params.amr_time_hours = 24

process importAccessions {
    label 'moshpit'

    cpus params.import_threads as int
    memory "${params.import_memory_gb} GB"
    time "${params.import_time_hours} h"

    input:
        path accessions_tsv

    output:
        path "accessions.qza", emit: accession_ids

    script:
        """
        qiime tools import \
            --type NCBIAccessionIDs \
            --input-path "${accessions_tsv}" \
            --output-path accessions.qza
        """
}

process fondueDownload {
    label 'moshpit'

    cpus params.fondue_threads as int
    memory "${params.fondue_memory_gb} GB"
    time "${params.fondue_time_hours} h"
    // q2-fondue runs fasterq-dump in /tmp. On Google Batch, specifying a
    // disk type mounts this dedicated scratch disk at /tmp.
    disk 200.GB, type: 'pd-standard'

    input:
        path accession_ids
        val fondue_email

    output:
        path "metadata.qza", emit: metadata
        path "single_reads.qza", emit: single_reads
        path "paired_reads.qza", emit: paired_reads
        path "failed_runs.qza", emit: failed_runs

    script:
        """
        qiime fondue get-all \
            --i-accession-ids ${accession_ids} \
            --p-email "${fondue_email}" \
            --p-threads ${task.cpus} \
            --o-metadata metadata.qza \
            --o-single-reads single_reads.qza \
            --o-paired-reads paired_reads.qza \
            --o-failed-runs failed_runs.qza \
            --verbose
        """

}

process assembleMegahit {
    label 'moshpit'

    cpus params.megahit_threads as int
    memory "${params.megahit_memory_gb} GB"
    time "${params.megahit_time_hours} h"
    // q2-assembly and MEGAHIT keep their working data under /tmp. Google
    // local SSDs are provisioned in 375 GB units and are suited to this
    // high-I/O assembly step.
    disk 375.GB, type: 'local-ssd'

    input:
        path reads

    output:
        path "contigs.qza", emit: contigs

    script:
    """
    qiime assembly assemble-megahit \
        --i-reads "${reads}" \
        --p-num-cpu-threads ${task.cpus} \
        --o-contigs contigs.qza \
        --verbose
    """
}

process predictGenesProdigal {
    label 'moshpit'

    cpus params.prodigal_threads as int
    memory "${params.prodigal_memory_gb} GB"
    time "${params.prodigal_time_hours} h"

    // q2-annotate unpacks the contigs artifact and writes multiple large
    // outputs under /tmp. Google local SSDs are provisioned in 375 GB units.
    disk 375.GB, type: 'local-ssd'

    input:
        path contigs

    output:
        path "loci.qza", emit: loci
        path "genes.qza", emit: genes
        path "proteins.qza", emit: proteins

    script:
        """
        qiime annotate predict-genes-prodigal \
          --i-seqs "${contigs}" \
          --p-mode meta \
          --o-loci loci.qza \
          --o-genes genes.qza \
          --o-proteins proteins.qza \
          --verbose
        """
}

process downloadAMRDB {
    label 'pathogenome'

    cpus params.amr_db_threads as int
    memory "${params.amr_db_memory_gb} GB"
    time "${params.amr_db_time_hours} h"

    output:
        path "amrfinderplus-db.qza", emit: amrfinderplus_db

    script:
        """
        qiime amrfinderplus fetch-amrfinderplus-db \
          --o-amrfinderplus-db amrfinderplus-db.qza \
          --verbose
        """
}

process AMRAnnotate {
    label 'pathogenome'

    cpus params.amr_threads as int
    memory "${params.amr_memory_gb} GB"
    time "${params.amr_time_hours} h"

    input:
        path amrfinderplus_db
        path contigs
        path proteins
        path loci

    output:
        path "amrfinderplus_annotations.qza", emit: annotations
        path "amrfinderplus_all_mutations.qza", emit: all_mutations
        path "amrfinderplus_genes.qza", emit: amr_genes
        path "amrfinderplus_proteins.qza", emit: amr_proteins

    script:
        """
        qiime amrfinderplus annotate \
          --i-amrfinderplus-db "${amrfinderplus_db}" \
          --i-sequences "${contigs}" \
          --i-proteins "${proteins}" \
          --i-loci "${loci}" \
          --p-threads ${task.cpus} \
          --o-amr-annotations amrfinderplus_annotations.qza \
          --o-amr-all-mutations amrfinderplus_all_mutations.qza \
          --o-amr-genes amrfinderplus_genes.qza \
          --o-amr-proteins amrfinderplus_proteins.qza \
          --verbose
        """
}


workflow {
    if (!params.fondue_email) {
        error "Missing required parameter: --fondue_email"
    }

    if (!params.accessions) {
        error "Missing required parameter: --accessions"
    }

    accessions_ch = channel.fromPath(params.accessions, checkIfExists: true)

    importAccessions(accessions_ch)

    fondueDownload(importAccessions.out.accession_ids, params.fondue_email)

    reads_ch = fondueDownload.out.paired_reads

    assembleMegahit(reads_ch)

    predictGenesProdigal(assembleMegahit.out.contigs)

    downloadAMRDB()

    AMRAnnotate(
        downloadAMRDB.out.amrfinderplus_db,
        assembleMegahit.out.contigs,
        predictGenesProdigal.out.proteins,
        predictGenesProdigal.out.loci
    )
}
