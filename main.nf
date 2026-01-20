nextflow.enable.dsl=2

params.input    = params.input  ?: "data/*_{R1,R2}.fastq.gz"
params.outdir   = params.outdir ?: "results"

process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"

    input:
        tuple val(sample_id), path(reads)

    output:
        buple val(sample_id), path("*_fastqc.zip"), path("*_fastqc.html")

    script:
    """
    fastqc -t ${task.cpus} ${reads}
    """
}

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    container "quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0"

    input:
        path(fastqpc_reports)

    output:
        path("multiqc_report.html")
        path("multiqc_data")

    script:
    """
    multiqc .
    """
}

workflow {
    // 1) Collect input reads
    Channel
        .fromPath(params.input, checkIfExists: true)
        .map { file ->
            // Extract sample_id from filename like SAMPLE_R1.fastq.gz or SAMPLE_R2.fastq.gz
            def name - file.getBaseName() // e.g. "SAMPLE_R1.fastq"
            def sample_id = name.replaceAll(/_R[12].*/, "")
            tuple(sample_id, file)
        }
        .groupTuple()
        .map { sample_id, files ->
            // files should be [R1, R2] (order not guaranteed)
            tuple(sample_id, files.sort { it.name })
        }
        .set { read_pairs_ch }

    // 2) FastQC on both reads for each sample
    fastqc_out = FASTQC(read_paris_ch)

    // 3) MultiQC over all FastQC outputs
    // FASTQC outputs tuples; MultiQC wants a folder of files.
    fastqc_files = fastqc_out
        .flatMap { sample_id, zip, html -> [zip, html] }
        .collect()

    MULTIQC(fastqc_files)
}