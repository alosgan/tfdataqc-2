/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CAT_FASTQ } from '../../../modules/nf-core/cat/fastq/main'
include { FASTQC as FASTQC_RAW } from '../../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIMMED } from '../../../modules/nf-core/fastqc/main'
include { FASTP } from '../../../modules/nf-core/fastp/main'
include { FASTP_LONG } from '../../../modules/local/fastplong/main'


workflow FASTQC_FASTP {
    take:
    ch_reads // channel: [ val(meta), [reads] ]
    fastp_save_trimmed_fail // boolean: true/false
    fastp_use_fastplong // boolean: true/false

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Make channels for concatenating fastq files
    ch_reads
        .branch { meta, fastqs ->
            single: fastqs.size() == 1
            return [meta, fastqs.flatten()]
            multiple: fastqs.size() > 1
            return [meta, fastqs.flatten()]
        }
        .set { ch_fastq }


    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ(
        ch_fastq.multiple
    ).reads.mix(ch_fastq.single).set { ch_filtered_reads }

    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

    //
    // MODULE: Run FASTQC_RAW
    //
    FASTQC_RAW(
        ch_filtered_reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect { it[1] })
    ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())

    // Change
    if (fastp_use_fastplong) {
        //
        // MODULE: Run fastplong
        //
        FASTP_LONG(
            ch_filtered_reads,
            [],
            false,
            fastp_save_trimmed_fail,
            false,
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTP_LONG.out.json.collect { it[1] })
        // Gather a list of the json files to throw into multiqc 
        ch_versions = ch_versions.mix(FASTP_LONG.out.versions.first())
        //
        // MODULE: Run FASTQC_TRIMMED
        //
        FASTQC_TRIMMED(
            FASTP_LONG.out.reads
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMMED.out.zip.collect { it[1] })
        ch_versions = ch_versions.mix(FASTQC_TRIMMED.out.versions.first())
    }
    else {
        //
        // MODULE: Run fastp
        //
        FASTP(
            ch_filtered_reads,
            [],
            false,
            fastp_save_trimmed_fail,
            false,
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect { it[1] })
        // Gather a list of the json files to throw into multiqc 
        ch_versions = ch_versions.mix(FASTP.out.versions.first())
        //
        // MODULE: Run FASTQC_TRIMMED
        //
        FASTQC_TRIMMED(
            FASTP.out.reads
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMMED.out.zip.collect { it[1] })
        ch_versions = ch_versions.mix(FASTQC_TRIMMED.out.versions.first())
    }

    if (fastp_use_fastplong) {
        reads = FASTP_LONG.out.reads
    }
    else {
        reads = FASTP.out.reads
    }

    emit:
    reads = reads
    multiqc_files = ch_multiqc_files
    version = ch_versions
}
