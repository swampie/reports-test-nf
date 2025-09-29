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
    path(pdb_file)
    
    output:
    path(pdb_file)
    


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
        Channel.of("${projectDir}/resources/AF-Q5VSL9-F1-model_v4.pdb")
    )
}
