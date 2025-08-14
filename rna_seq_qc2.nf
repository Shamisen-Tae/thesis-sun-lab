#!/usr/bin/env nextflow

params.reads = "${params.fasta_dir}/*.fastq.gz"
params.outdir = "${params.arc_dir}/results"
params.pairs = "${params.fasta_dir}/*_R{1,2}_001.fastq.gz"
params.genomeDir = "${params.arc_dir}/hg38_index"
params.trimmed = "${params.outdir}/fastp/*_R{1,2}.trimmed.fastq.gz"

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
        --thread ${task.cpus} 
    """
    // fastp use --thread while fastqc use --threads :((
    // fastp: downgrade to v0.20.1 to in conda env!!
}


// align reads with STAR
// mark duplicate reads with Picard from STAR 
// assign reads to genes with featureCounts with bam file from Picard
process star_align {
    publishDir "${params.outdir}/star_alignment", mode: 'move'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id),
        path("${sample_id}.Aligned.sortedByCoord.out.bam"),
        path("${sample_id}.Log.final.out")

    script:
    """
    STAR --runThreadN ${task.cpus} \
        --genomeDir ${params.genomeDir} \
        --readFilesIn ${reads.join(' ')} \
        --readFilesCommand zcat \
        --outFileNamePrefix ${sample_id}. \
        --outSAMtype BAM SortedByCoordinate \
        --outReadsUnmapped Fastx \
        --outFilterMatchNminOverLread 0.50 \
        --outFilterScoreMinOverLread 0.50 \
        --outSAMmapqUnique 40 \
        --outFilterMultimapNmax 1
    """
}

process mark_duplicates {
    publishDir "${params.outdir}/deduplicate_bam", mode: 'move'

    input:
    path bam

    output:
    tuple val(bam.baseName), path("${bam.baseName}.dedup.bam"), path("${bam.baseName}.dedup.metrics.txt")

    script:
    """
    picard MarkDuplicates \
        I=${bam} \
        O=${bam.baseName}.dedup.bam \
        M=${bam.baseName}.dedup.metrics.txt \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=SILENT
    """
}

process feature_counts {
    publishDir "${params.outdir}/featurecounts", mode: 'move'

    input:
    path dedup_bam

    output:
    path "${dedup_bam.baseName}.counts.txt"

    script:
    """
    featureCounts -T ${task.cpus} \
        -p \
        -t exon \
        -g gene_id \
        -a "${params.genomeDir}/Homo_sapiens.GRCh38.110.gtf" \
        -o ${dedup_bam.baseName}.counts.txt \
        ${dedup_bam}
    """
}



workflow {
    // read_fastas = Channel
    //     .fromPath(params.reads).collate(params.batch_size)
    // read_fastps = Channel.fromFilePairs(params.pairs, flat: true) 
    //     .map { sample_id, r1, r2 -> tuple(sample_id, [r1, r2]) }
    
    //fastqc(read_fastas)
    //fastp(read_fastps)

    // read_trimmed = Channel.fromFilePairs(params.trimmed, flat: true)
  	//    .map { sample_id, r1, r2 -> tuple(sample_id, [r1, r2]) }   
    

   //params.aligned = "${params.arc_dir}/results/star_alignment/*.bam"
   //read_aligned = Channel.fromPath(params.aligned)

   //mark_duplicates(read_aligned)
   params.deduped = "${params.outdir}/deduplicate_bam/*.bam"

   deduped = Channel.fromPath(params.deduped)

   feature_counts(deduped)


}
