# Usage: snakemake --snakefile getToolVersions.smk --printshellcmds --jobs 1 --latency-wait 60 --use-singularity --singularity-args "--bind $(pwd)"
# NB: GeCKO.sif must already be available in GeCKO/utils/singularity_image

path_to_snakefile = workflow.snakefile
snakefile_dir = path_to_snakefile.rsplit('/', 1)[0]
GeCKO_image = os.path.abspath(os.path.join(snakefile_dir, "../singularity_image/GeCKO.sif"))



rule All:
    input:
        "tool_versions.txt",
        "python_lib_versions.txt"


rule check_tools:
    input:
        "CONFIG/tools.txt"
    output:
        "tool_versions.txt"
    singularity:
        GeCKO_image
    shell:
        """
        SCRIPTS/getToolVersions.sh {input} {output}
        """

rule check_python_libraries:
    input:
        "CONFIG/python_libs.txt"
    output:
        "python_lib_versions.txt"
    singularity:
        GeCKO_image
    shell:
        "python SCRIPTS/getLibVersions.py -i {input} -o {output}"
