# Pathogenome AMR workflow

**This workflow is currently a beta release.**

This Nextflow workflow starts with a list of NCBI or ENA sequencing-run accessions and looks for antimicrobial resistance (AMR) genes in their assembled genomes.
It is configured to run on the Seqera Platform and uses the [pathogenome](https://library.qiime2.org/quickstart/pathogenome) distribution of QIIME 2 under the hood.

## What the workflow does

1. Imports accession IDs.
2. Downloads the sequencing data with `q2-fondue`.
3. Assembles paired-end reads into contigs with MEGAHIT.
4. Predicts genes and proteins with Prodigal.
5. Downloads the AMRFinderPlus database and annotates the predicted genes and proteins for AMR-related features.

Only paired-end reads are passed to the assembly step.

## Required inputs

### Accession file

Provide a tab-separated text file with a `sample-id` header and one NCBI or ENA run accession per row:

```text
sample-id
ERR00000000
ERR00000001
```

Upload this file to a location available in Seqera Data Explorer, or add it as a Seqera dataset.

### Contact email

The `q2-fondue` step requires an email address (because NCBI requires one).
Supply one on the launch page on Seqera.

## Launching on Seqera

1. Open **Launchpad** in your Seqera workspace and add or select this pipeline.
2. Choose the `v0.1.0-beta` revision.
3. Select your Google Batch compute environment and Google Cloud Storage work directory.
4. Leave **Config profile** empty and set **Pipeline schema** to **Repository default**.
5. In **Run parameters**, select the accession TSV for **Accessions** and enter your address for **Fondue email**.
6. Optionally adjust the Fondue and MEGAHIT thread counts.
7. Select **Launch**.

## Finding the results

After the run finishes, open **Tasks**, select the `AMRAnnotate` task, and open the **Data Explorer** tab.
The QIIME 2 outputs are:

- `amrfinderplus_annotations.qza`
- `amrfinderplus_all_mutations.qza`
- `amrfinderplus_genes.qza`
- `amrfinderplus_proteins.qza`

Intermediate artifacts such as the reads downloaded with `q2-fondue`, or the contigs from assembly can be found in the Data Explorer tabs of the corresponding tasks.
