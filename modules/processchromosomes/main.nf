//Processing chromosomes using python script
//Python script checks for polyA sites near end of transcripts that have FANTOM support 
process PROCESSCHROMOSOMES {
    //publishDir 'outputs/processedChrs', mode: 'copy', overwrite: true
    input:
    tuple path(fantom), path(encode), path(polyA), path(human),  val(id), val(chr)

    output:
    tuple val(id), val(chr), path('output_*')
    

    """
    polyAFantom3primeChecker.py ${human} ${fantom} ${encode} ${polyA} output ${chr}
    """
}