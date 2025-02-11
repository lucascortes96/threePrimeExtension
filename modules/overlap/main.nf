process OVERLAP {
    publishDir 'outputs/overlap', mode: 'copy', overwrite: true
    input:
    path result
    path overlap 
    path gff
    output:
    path 'overlapped.csv'

    """
    python3 ${overlap} ${result} ${gff} overlapped.csv
    """
}