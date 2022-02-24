reference = config["REFERENCE"]


rule index:
    input:
        reference
    output:
        reference+".bwt.2bit.64"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "bwa-mem2 index {input}"

