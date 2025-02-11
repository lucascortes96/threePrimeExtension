//Returns split transcript, polyA and chr for next process
process SPLITCHROMOSOMES {
    input:
    tuple val(id), path(polyA), path(fantom), path(encode), path(human), val (chr)
    
    output:
    tuple path ('split_fantom_*'), path('split_encode_*'), path('split_polyA_*'), path('split_human_*'),  val(id), val(chr)
    //val (chr), emit:chr
    

    """
    awk -v chr="${chr}" 'NR==1 {header=\$0; next} \$1 == chr {if (!printed) {print header > "split_fantom_" chr ".txt"; printed=1} print >> "split_fantom_" chr ".txt"}' ${fantom}
    awk -v chr="${chr}" 'NR==1 {header=\$0; next} \$1 == chr {if (!printed) {print header > "split_encode_" chr ".txt"; printed=1} print >> "split_encode_" chr ".txt"}' ${encode}
    awk -v chr="${chr}" '
    BEGIN { FS = OFS = "\t" }
    NR==1 {
        header = "Chromosome\tStart\tEnd\tID\tScore\tStrand\tSignal\tNumber\tMean RPM\tType\tMotifs"
        next
    }
    \$1 == chr {
        if (!printed) {
            print header > "split_polyA_" chr ".txt"
            printed = 1
        }
        print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11 >> "split_polyA_" chr ".txt"
    }
    ' ${polyA}
    awk -v chr="${chr}" 'NR==1 {header=\$0; next} \$1 == chr {if (!printed) {print header > "split_human_" chr ".txt"; printed=1} print >> "split_human_" chr ".txt"}' ${human}
    """
}
