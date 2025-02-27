process CLEANUP {
    publishDir 'outputs/cleaned', mode: 'copy', overwrite: true
    input:
    tuple val(id), path(result)
    output:
    path "${id}_extended_transcripts_final.csv", emit: csv
    val (id), emit:id
    path "${id}_extended_transcripts_threePrimeExtendStats.txt"
    path "${id}_extended_transcripts_threePrimeExtendPlot.png"

    """
    finalFilterandStats.py ${result} "${id}_extended_transcripts.csv"
    """
}