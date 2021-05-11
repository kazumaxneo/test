#!/usr/bin/env nextflow

Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs_ch }

/*
 * Step 1. Adaptor trimming and qulaity trimming
 */
process preprocessing {

    tag "$pair_id"
    publishDir params.outdir1, mode: 'copy'
    input:
    tuple val(pair_id), path(reads) from read_pairs_ch
    output:
    set pair_id, "${pair_id}_1_val_1.fq.gz" into read1_ch
    set pair_id, "${pair_id}_2_val_2.fq.gz" into read2_ch
    set pair_id, "*fastqc.html"

    """
    trim_galore --paired ${pair_id}_1.fastq.gz ${pair_id}_2.fastq.gz --fastqc --length 20 -q 30
    """
}

/*
 * Step 2. BwaIndexing
 */
process bwaindex {
    tag "$genome.baseName"
    input:
    path genome from params.genome
    output:
    path "bwa*" into bwaindex_ch

    """
    bwa index -a is -p bwaindex ${genome} 
    """
}

/*
 * Step 3. Mapping
 */
process mapping {
    tag "$pair_id"
    publishDir params.outdir2, mode: 'copy'

    input:
    path index from bwaindex_ch
    tuple val(pair_id), path(reads1) from read1_ch
    tuple val(pair_id), path(reads2) from read2_ch

    output:
    set pair_id, "${pair_id}.bam" into bam_ch

    """
    bwa mem -t ${task.cpus} bwaindex ${pair_id}_1_val_1.fq.gz ${pair_id}_2_val_2.fq.gz \
    |samtools sort -O BAM -@ ${task.cpus} - > ${pair_id}.bam
    """
}

/*
 * Step X. FastaIndexing
 */
process referenceindex {
    tag "$genome.baseName"
    input:
    path genome from params.genome
    output:
    path "${genome}" into fasta_ch
    path "${genome}.fai" into fastafai_ch

    """
    samtools faidx ${genome} 
    """
}

/*
 * Step X. variant calling
 */
process variantcall {
    tag "$pair_id"
    publishDir params.outdir3, mode: 'copy'

    input:
    path fastaindex from fasta_ch
    path fastafaiindex from fastafai_ch
    tuple val(pair_id), path(bam) from bam_ch

    output:
    set pair_id, "${pair_id}_freebayes.vcf"

    """
    freebayes -u -f ${fastaindex} $bam > ${pair_id}_freebayes.vcf
    """
}
