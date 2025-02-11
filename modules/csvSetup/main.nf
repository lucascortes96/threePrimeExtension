process CSVSETUP {
    input:
    val ch_input

    output:
    path 'updated_output.csv'
    

    script:
    """
    echo "Received input: '${ch_input}'"
    csvSetup.py '${ch_input.join(',')}'
    """
}