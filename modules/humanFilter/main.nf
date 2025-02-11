process HUMANFILTER{
    input:
    tuple val(id), path(human), path(polyA), path(fantom), path(encode)
    path readthroughs
    val single_exon
    output:
    tuple val(id), path(polyA), path(fantom), path(encode), path ('noReadthroughProteinCoding.gff3')
    """
    humanFilter.py ${human} 'noReadthroughProteinCoding.gff3' ${readthroughs} ${single_exon}
    """
}