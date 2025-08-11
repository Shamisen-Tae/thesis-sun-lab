#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = "${params.fasta_dir}/*.fastq.gz"
params.outdir = "${params.arcdir}/results"
params.pairs = "data/*_R{1,2}.fastq.gz"


process fastqc {
    publishDir "${params.outdir}/fastqc", mode: 'move'

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


process fastp {
    tag "$sample_id"
    publishDir "${params.outdir}/fastp", mode: 'move'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_R1.trimmed.fastq.gz"), path("${sample_id}_R2.trimmed.fastq.gz"), path("${sample_id}_fastp.html"), path("${sample_id}_fastp.json")

    script:
    """
    echo "Running fastp on: ${reads.join(' ')}"
    fastp \
        --in1 "${reads[0]}" \
        --in2 "${reads[1]}" \
        --out1 "${sample_id}_R1.trimmed.fastq.gz" \
        --out2 "${sample_id}_R2.trimmed.fastq.gz" \
        --length_required 25 \
        --cut_tail \
        --cut_tail_window_size 4 \
        --cut_tail_mean_quality 20 \
        --disable_quality_filtering \
        --overrepresentation_analysis \
        --html "${sample_id}_fastp.html" \
        --json "${sample_id}_fastp.json" \
        --threads ${task.cpus}
    """
}


workflow {
    read_fastas = Channel
        .fromPath(params.reads).collate(params.batch_size)
    read_fastps = Channel.fromFilePairs(params.pairs, flat: true)
      
    //fastqc(read_fastas)
    fastp(read_fastps)

    

}