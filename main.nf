#!/usr/bin/env nextflow

params.accessions = null
params.fondue_email = null

params.fondue_threads = 8
params.megahit_threads = 12


process importAccessions {
    label 'moshpit'

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

    accessions_path = params.accessions ?: "${projectDir}/assets/accessions.tsv"
    accessions_ch = channel.fromPath(accessions_path, checkIfExists: true)

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
