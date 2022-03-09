rule all:
    input:
    output:
        res.txt
    params:
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "gatk --help > res.txt"
