/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This borrows code heavily from nf-core/taxprofiler.  Though is MUCH simpler in ambition.
*/
include { UNTAR } from '../../modules/nf-core/untar/main'
include { KRAKEN2_KRAKEN2 } from '../../modules/nf-core/kraken2/kraken2/main'
include { KRAKENTOOLS_KREPORT2KRONA } from '../../modules/nf-core/krakentools/kreport2krona/main'
include { KRONA_CLEANUP } from '../../modules/local/krona_cleanup'
include { KRONA_KTIMPORTTEXT } from '../../modules/nf-core/krona/ktimporttext/main' 

workflow KRAKEN2_WF{
    take:
    kraken2_db             // string: path to kraken2db either tar.gz or folder
    ch_trimmed_reads               // channel: [ val(meta), [reads] ]

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_krona_text = Channel.empty()

    ch_kraken_db = Channel.fromPath(kraken2_db)

    ch_kraken_db | map{ k_path ->
        [[id:k_path.getSimpleName()],file(k_path)]
    }
    | branch{ meta, k_path ->
        untar: k_path.name.endsWith(".tar.gz")
        skip:true
    } | set { ch_kraken_db }
    
    //
    // UNTAR if required
    //
    UNTAR( ch_kraken_db.untar )
    ch_versions = ch_versions.mix( UNTAR.out.versions.first() )
    ch_final_dbs = ch_kraken_db.skip.mix(UNTAR.out.untar).map{ meta, db ->
        [db]
    }
    
    // Prepare channel for kraken2 so we've got same number of DBs and reads
    ch_trimmed_reads.combine(ch_final_dbs).multiMap{meta, read_path, database_path ->
        reads:[meta, read_path]
        db:[database_path]
    } | set { ch_trimmed_reads }
    
    //
    //  KRAKEN2 analysis of reads
    //
    KRAKEN2_KRAKEN2(
        ch_trimmed_reads.reads,
        ch_trimmed_reads.db,
        params.kraken2_save_output_fastqs,
        params.kraken2_save_reads_assignment,
    )
    ch_multiqc_files       = ch_multiqc_files.mix( KRAKEN2_KRAKEN2.out.report.collect{it[1]} )
    ch_versions            = ch_versions.mix( KRAKEN2_KRAKEN2.out.versions.first() )

    //
    //  Convert from KRAKEN2 report to KRONA formatted
    //
    KRAKENTOOLS_KREPORT2KRONA ( KRAKEN2_KRAKEN2.out.report )
    ch_krona_text = ch_krona_text.mix( KRAKENTOOLS_KREPORT2KRONA.out.txt )
    ch_versions = ch_versions.mix( KRAKENTOOLS_KREPORT2KRONA.out.versions.first() )

    //
    // Remove taxonomy level annotations from the Krona text files
    //
    KRONA_CLEANUP( ch_krona_text )
    ch_cleaned_krona_text = KRONA_CLEANUP.out.txt
    ch_versions = ch_versions.mix( KRONA_CLEANUP.out.versions.first() )

    // Combine into single list so all are in on html file
    ch_krona_text_for_import = ch_cleaned_krona_text
        .map{
            [[id:'kraken2'], it[1]]}
        .groupTuple()
        .map{
            meta, files ->
            [meta, files.sort { it.name.toLowerCase() }]
        }

    //
    // Convert Krona text files into html Krona visualizations
    //
    KRONA_KTIMPORTTEXT( ch_krona_text_for_import )
    ch_krona_html = KRONA_KTIMPORTTEXT.out.html
    ch_versions = ch_versions.mix( KRONA_KTIMPORTTEXT.out.versions.first() )
    
    emit:
    multiqc_files = ch_multiqc_files
    version = ch_versions
    kraken_report = KRAKEN2_KRAKEN2.out.report
    kraken_classified_reads_assignment = KRAKEN2_KRAKEN2.out.classified_reads_assignment
}