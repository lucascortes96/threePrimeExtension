process GENERALEXONMATCHER{
    memory '40 GB'
    cpus 2
    input:
    tuple val(id), path(polyA), path(fantom), path(encode), path (human)
    val single_exon
    output:
    tuple val(id), path(polyA), path("matched_fantom_blocks.gff"), path("matched_encode_blocks.gff"), path("filtered_matched_human_exons.gff")
    """
    generalExonMatcher.py ${human} ${fantom} ${encode} ${single_exon}  .
    """
}