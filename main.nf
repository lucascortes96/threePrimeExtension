include { validateParameters; paramsSummaryLog; paramsSummaryMap; samplesheetToList } from 'plugin/nf-schema'

// READING IN THE CHROMOSOMES
chromosomes = Channel.from('1', '2',
'3', '4', '5', '6', '7', '8', '9', 
'10','11', '12', '13', '14', '15', 
'16', '17', '18', '19', '20','21', 
'22', 'X', 'Y') // ADD X AND Y BACK IN LATER

// Import modules 
include { SPLITCHROMOSOMES } from './modules/splitchromosomes'
include { PROCESSCHROMOSOMES } from './modules/processchromosomes'
include { CATALL } from './modules/catall'
include { CLEANUP } from './modules/cleanup'
include { OVERLAP } from './modules/overlap'
include {THREEPRIMEGRAB} from './modules/threePrimeGrab'
include {GENERALEXONMATCHER} from './modules/generalExonMatcher'
include {HUMANFILTER} from './modules/humanFilter'
include {CSVSETUP} from './modules/csvSetup'


// Validate input parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)

workflow {
    //READING IN PARAMS
    outputDir = file(params.outputDir)
    human = file(params.human)
    polyA = file(params.polyA)
    fantom = file(params.fantom)
    encode = file(params.encode)
    cleanup = file(params.cleanup)
    referenceGff = file(params.referenceGff)
    readThroughs = file(params.readThroughs)
    
    
    // Assuming 'input' is the parameter name for your CSV file in your nextflow_schema.json
    ch_input = Channel.fromList(samplesheetToList(params.input_csv, "${projectDir}/assets/input_schema.json")).collect().view()
    
    //.collect().view()
    generalChannel = CSVSETUP(ch_input)
    generalChannel.view()
    generalChannel = generalChannel.splitCsv()
    generalChannel.view()
    //RUNNNING WORKFLOW 
    humanOut = HUMANFILTER( generalChannel, readThroughs, params.single_exon)
    humanOut.collect().view()

    threePrimeOut = THREEPRIMEGRAB(humanOut).view()
    generalOut = GENERALEXONMATCHER(threePrimeOut, params.single_exon)
    generalChromosomes = generalOut.combine(chromosomes).view()
    splitChrs = SPLITCHROMOSOMES(generalChromosomes) 
    processChrOut = PROCESSCHROMOSOMES(splitChrs)
    .groupTuple(by: 0)  // Group by ID (index 0)
    .map { id, chrs, files -> tuple(id, files.flatten()) }  // Flatten the list of files
    .view()

    catted = CATALL(processChrOut)
    CLEANUP(catted)  
    //OVERLAP(CLEANUP.out, overlap, referenceGff)
}