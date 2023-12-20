nextflow.enable.dsl = 2

params.outdir = "results"
params.delay = 5

process MULTIQC {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(multiqc_html) 
    val(step)

    output:
    path("step_*/*.html")

    script:
    """
    echo "Copying MultiQC reports"
    mkdir step_$step
    cp $multiqc_html step_$step/
    """
}

process REPORTS {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val multiqc_reports
    path(bin_depths_summary_tsv)
    path(bin_summary_tsv)
    path(CAPES_S7_log)
    path(execution_trace_txt)
    path(kraken2_report_txt)
    path(pipeline_dag_svg)
    path(binDepths_png)
    path(krona_html)
    path(transposed_tex)
    path(transposed_tsv)
    path(transposed_txt)
    path(genome_fasta)
    path(all_sites_fas)
    path(baits_bed)
    path(genome_dict)
    path(genome_fasta_fai)
    path(genome_gff3)
    path(genome_gtf)
    path(genome_sizes)
    path(proteome_fasta)
    path(test_baserecalibrator_table)
    path(test_bed)
    path(test_bedgraph)
    path(test_bigwig)
    path(test_paired_end_bam)
    path(test_single_end_bam_readlist_txt)
    path(test_vcf)
    path(test_vcf_gz_tbi)
    path(test2_targets_tsv_gz)
    path(transcriptome_paf)
    path(report_pdf)
    path(samplesheet_csv)
    path(report_json)


    output:
    tuple path(bin_depths_summary_tsv),
        path(bin_summary_tsv),
        path(CAPES_S7_log),
        path(execution_trace_txt),
        path(kraken2_report_txt),
        path(pipeline_dag_svg),
        path(binDepths_png),
        path(krona_html),
        path(transposed_tex),
        path(transposed_tsv),
        path(transposed_txt),
        path(genome_fasta),
        path(report_pdf),
        path(all_sites_fas),
        path(baits_bed),
        path(genome_dict),
        path(samplesheet_csv),
        path(genome_fasta_fai),
        path(genome_gff3),
        path(genome_gtf),
        path(genome_sizes),
        path(proteome_fasta),
        path(test_baserecalibrator_table),
        path(test_bed),
        path(test_bedgraph),
        path(test_bigwig),
        path(test_paired_end_bam),
        path(test_single_end_bam_readlist_txt),
        path(test_vcf),
        path(test_vcf_gz_tbi),
        path(test2_targets_tsv_gz),
        path(transcriptome_paf),
        path(report_json)



    script:
    """
    echo "Sleeping ${params.delay} seconds"
    sleep ${params.delay}
    echo "Copying all resource files to results directory for testing!"
    """

}

workflow {
    MULTIQC(Channel.value("${projectDir}/resources/MultiQC Report.html"), Channel.from(1..4))
    
    REPORTS(
        MULTIQC.out,
        Channel.of("${projectDir}/resources/bin_depths_summary.tsv"),
        Channel.of("${projectDir}/resources/bin_summary.tsv"),
        Channel.of("${projectDir}/resources/CAPES_S7.log"),
        Channel.of("${projectDir}/resources/execution_trace_2021-07-29_07-12-59.txt"),
        Channel.of("${projectDir}/resources/kraken2_report.txt"),
        Channel.of("${projectDir}/resources/pipeline_dag_2021-07-29_07-12-59.svg"),
        Channel.of("${projectDir}/resources/SPAdesHybrid-CAPES_S7-binDepths.heatmap.png"),
        Channel.of("${projectDir}/resources/taxonomy.krona.html"),
        Channel.of("${projectDir}/resources/transposed_report.tex"),
        Channel.of("${projectDir}/resources/transposed_report.tsv"),
        Channel.of("${projectDir}/resources/transposed_report.txt"),
        Channel.of("${projectDir}/resources/genome.fasta"),
        Channel.of("${projectDir}/resources/all_sites.fas"),
        Channel.of("${projectDir}/resources/baits.bed"),
        Channel.of("${projectDir}/resources/genome.dict"),
        Channel.of("${projectDir}/resources/genome.fasta.fai"),
        Channel.of("${projectDir}/resources/genome.gff3"),
        Channel.of("${projectDir}/resources/genome.gtf"),
        Channel.of("${projectDir}/resources/genome.sizes"),
        Channel.of("${projectDir}/resources/proteome.fasta"),
        Channel.of("${projectDir}/resources/test.baserecalibrator.table"),
        Channel.of("${projectDir}/resources/test.bed"),
        Channel.of("${projectDir}/resources/test.bedgraph"),
        Channel.of("${projectDir}/resources/test.bigwig"),
        Channel.of("${projectDir}/resources/test.paired_end.bam"),
        Channel.of("${projectDir}/resources/test.single_end.bam.readlist.txt"),
        Channel.of("${projectDir}/resources/test.vcf"),
        Channel.of("${projectDir}/resources/test.vcf.gz.tbi"),
        Channel.of("${projectDir}/resources/test2.targets.tsv.gz"),
        Channel.of("${projectDir}/resources/transcriptome.paf"),
        Channel.of("${projectDir}/resources/report.pdf"),
        Channel.of("${projectDir}/resources/nfcore_chipseq110_samplesheet_test_full_6cols.csv"),
        Channel.of("${projectDir}/resources/report.json")
    )
}
