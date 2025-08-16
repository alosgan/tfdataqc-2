process FASTP_LONG {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/fastplong:0.2.2--heae3180_0'
        : 'biocontainers/fastplong:0.2.2--heae3180_0'}"

    input:
    tuple val(meta), path(reads)
    path adapter_fasta
    // File in FASTA format containing possible adapters to remove.
    val discard_trimmed_pass
    // Specify true to not write any reads that pass trimming thresholds. | This can be used to use fastp for the output report only.
    val save_trimmed_fail
    // Specify true to save files that failed to pass trimming thresholds ending in \*.fail.fastq.gz
    val save_merged

    output:
    tuple val(meta), path('*.fastp.fastq.gz'), optional: true, emit: reads
    tuple val(meta), path('*.json'), emit: json
    tuple val(meta), path('*.html'), emit: html
    tuple val(meta), path('*.log'), emit: log
    tuple val(meta), path('*.fail.fastq.gz'), optional: true, emit: reads_fail
    tuple val(meta), path('*.merged.fastq.gz'), optional: true, emit: reads_merged
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def adapter_list = adapter_fasta ? "--adapter_fasta ${adapter_fasta}" : ""
    def fail_fastq = save_trimmed_fail ? "--failed_out ${prefix}.fail.fastq.gz" : ""
    def out_fq1 = discard_trimmed_pass ? "" : "--out ${prefix}.fastp.fastq.gz"


    // Only running in one mode currenly -- unless I encounter more complex cases
    """
    [ ! -f  ${prefix}.fastq.gz ] && ln -sf ${reads} ${prefix}.fastq.gz

    fastplong \\
        --in ${prefix}.fastq.gz \\
        ${out_fq1} \\
        --thread ${task.cpus} \\
        --json ${prefix}.fastp.json \\
        --html ${prefix}.fastp.html \\
        ${adapter_list} \\
        ${fail_fastq} \\
        ${args} \\
        2> >(tee ${prefix}.fastp.log >&2)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastplong: \$(fastplong --version 2>&1 | sed -e "s/fastplong //g")
    END_VERSIONS
    """
}