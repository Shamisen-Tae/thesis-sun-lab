#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = "${params.fasta_dir}/*.fastq.gz"
params.outdir = "${params.arcdir}/results/fastqc"


process fastqc {
    publishDir params.outdir, mode: 'move'

    input:
    path reads

    output:
    path "*.html"
    path "*.zip"

    tag { "fastqc_batch_${task.index}" }

    script: 
    """
    echo "Running FastQC on: ${reads.join(' ')}"
    fastqc --threads ${task.cpus} ${reads.join(' ')} --outdir .
    """
}

workflow {
    read_fastas = Channel
        .fromPath(params.reads).collate(params.batch_size)
        
    fastqc(read_fastas)
}