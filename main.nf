#!/usr/bin/env nextflow

params.accessions = null
params.outdir = "results"

params.fondue_threads = 8
params.megahit_threads = 12


process importAccessions {
    label 'moshpit'

    publishDir "${params.outdir}/accessions", mode: "copy"

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

    publishDir "${params.outdir}/fondue", mode: 'copy'

    input:
        path accession_ids

    output:
        path "metadata.qza", emit: metadata
        path "single_reads.qza", emit: single_reads
        path "paired_reads.qza", emit: paired_reads
        path "failed_runs.qza", emit: failed_runs

    script:
        """
        qiime fondue get-all \
            --i-accession-ids ${accession_ids} \
            --p-email "REMOVED_EMAIL" \
            --p-threads ${task.cpus} \
            --o-metadata metadata.qza \
            --o-single-reads single_reads.qza \
            --o-paired-reads paired_reads.qza \
            --o-failed-runs failed_runs.qza \
            --verbose
        """

}

process exportFondueDiagnostics {
    label 'moshpit'

    publishDir "${params.outdir}/diagnostics/fondue",
        mode: "copy",
        saveAs: { filename ->
            filename == "paired_reads.qza" ? null : filename
        }

    input:
        path metadata
        path failed_runs
        path paired_reads

    output:
        path "metadata", emit: metadata
        path "failed_runs", emit: failed_runs
        path "paired_reads_manifest.csv", emit: paired_reads_manifest
        path "paired_reads_archive_listing.txt", emit: paired_reads_archive_listing
        path paired_reads, emit: paired_reads

    script:
        """
        qiime tools export \
            --input-path "${metadata}" \
            --output-path metadata

        qiime tools export \
            --input-path "${failed_runs}" \
            --output-path failed_runs

        unzip -p "${paired_reads}" '*/data/MANIFEST' \
            > paired_reads_manifest.csv
        unzip -l "${paired_reads}" \
            > paired_reads_archive_listing.txt
        """
}

process assembleMegahit {
    label 'moshpit'

    cpus params.megahit_threads as int

    publishDir "${params.outdir}/assembly", mode: "copy"

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

process exportAssemblyDiagnostics {
    label 'moshpit'

    publishDir "${params.outdir}/diagnostics/assembly",
        mode: "copy",
        saveAs: { filename ->
            filename == "contigs.qza" ? null : filename
        }

    input:
        path contigs

    output:
        path contigs, emit: contig_artifact
        path "contigs", emit: exported_contigs
        path "contig_stats.tsv", emit: contig_stats

    script:
        """
        qiime tools export \
            --input-path "${contigs}" \
            --output-path contigs

        printf 'file\\tcontigs\\tbases\\n' > contig_stats.tsv
        find contigs -type f -name '*.fa' -exec \
            awk 'BEGIN { n=0; bp=0 } /^>/ { n++; next } { bp += length(\$0) } END { printf "%s\\t%d\\t%d\\n", FILENAME, n, bp }' {} \\; \
            >> contig_stats.tsv
        """
}

process predictGenesProdigal {
    label 'moshpit'

    publishDir "${params.outdir}/gene_prediction", mode: "copy"

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

    publishDir "${params.outdir}/amrfinderplus_db", mode: "copy"

    output:
        path "amrfinderplus-db.qza", emit: amrfinderplus_db

    script:
        """
        git clone https://github.com/bokulich-lab/q2-amrfinderplus q2-amrfinderplus
        cd q2-amrfinderplus
        git checkout main
        python -m pip install -e .
        qiime dev refresh-cache
        cd -

        qiime amrfinderplus fetch-amrfinderplus-db \
          --o-amrfinderplus-db amrfinderplus-db.qza \
          --verbose
        """
}

process AMRAnnotate {
    label 'pathogenome'

    publishDir "${params.outdir}/05_amrfinderplus", mode: "copy"

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
          --o-amr-annotations amrfinderplus_annotations.qza \
          --o-amr-all-mutations amrfinderplus_all_mutations.qza \
          --o-amr-genes amrfinderplus_genes.qza \
          --o-amr-proteins amrfinderplus_proteins.qza \
          --verbose
        """
}


workflow {
    accessions_path = params.accessions ?: "${projectDir}/assets/accessions.tsv"
    accessions_ch = channel.fromPath(accessions_path, checkIfExists: true)

    importAccessions(accessions_ch)

    fondueDownload(importAccessions.out.accession_ids)

    exportFondueDiagnostics(
        fondueDownload.out.metadata,
        fondueDownload.out.failed_runs,
        fondueDownload.out.paired_reads
    )

    reads_ch = exportFondueDiagnostics.out.paired_reads

    assembleMegahit(reads_ch)

    exportAssemblyDiagnostics(assembleMegahit.out.contigs)

    predictGenesProdigal(exportAssemblyDiagnostics.out.contig_artifact)

    downloadAMRDB()

    AMRAnnotate(
        downloadAMRDB.out.amrfinderplus_db,
        exportAssemblyDiagnostics.out.contig_artifact,
        predictGenesProdigal.out.proteins,
        predictGenesProdigal.out.loci
    )
}
