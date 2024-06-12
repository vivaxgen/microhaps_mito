# get the vivaxGEN ngs-pipeline base directory
NGS_PIPELINE_BASE = config['NGS_PIPELINE_BASE']

IDs = read_files.keys()

rule final:
    input:
        "result.txt"

rule something:
    output:
        "result.txt"
    shell:
        "echo 'Hello World!' > result.txt"
