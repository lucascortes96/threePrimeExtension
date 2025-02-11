process THREEPRIMEGRAB {
    input:
    tuple val(id), path(polyA), path(fantom), path(encode), path (noReadThrough) 
    output:
    tuple val(id), path(polyA), path(fantom), path(encode), path ("grabbedThreePrimehg38.gff")
    """
    threePrimeGrab.py ${noReadThrough} grabbedThreePrimehg38.gff
    """
}