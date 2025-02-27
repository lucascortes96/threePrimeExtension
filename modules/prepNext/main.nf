
process PREPNEXT {
    publishDir 'outputs/cleaned', mode: 'copy', overwrite: true
    input:
    path (csv)
    val (id)
    output:
    path "${id}_nextRun.gff"

    """
    prepNext.py ${csv} ${id}
    """
}