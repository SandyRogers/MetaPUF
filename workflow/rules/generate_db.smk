import os

import pandas as pd

ASSEMBLIES_FOLDER = os.path.join(config["outputdir"], "assemblies")
SAMPLEINFO_FILE    = config["sample_metadata"]
SAMPLEINFO_FILE_FINAL = os.path.join(ASSEMBLIES_FOLDER, "sample_info_final.csv")
sample_info     = pd.read_csv(SAMPLEINFO_FILE, sep=',')
Assemblies      = sample_info['Assembly'].drop_duplicates().to_list()
OUTPUT_FILE = expand("{af}/{aname}_contig_info.txt", af=ASSEMBLIES_FOLDER, aname=Assemblies)
PROTEIN_FILE = expand("{af}/{aname}.faa.gz", af=ASSEMBLIES_FOLDER, aname=Assemblies)
CLUSTER_REPORT = os.path.join(ASSEMBLIES_FOLDER, "databases/cluster_report.txt")
STUDY       = os.environ.get("STUDY", config["parameters"]["study"])
OUTPUTDIR   = config["outputdir"]

rule generate_db:
    input:
        sample_metadata=SAMPLEINFO_FILE
    output:
        contigs_dir=OUTPUT_FILE,
        protein_file=PROTEIN_FILE,
        sample_info_final=SAMPLEINFO_FILE_FINAL,
        cluster_rpt=CLUSTER_REPORT
    params:
        study=STUDY,
        input_dir=config["parameters"]["input_dir"],
        ver=config["parameters"]["mgnify_version"],
        output_dir=OUTPUTDIR,
        db_size=config["parameters"]["db_size"]
    log:
        expand("logs/{iname}_db_generate.log",iname=STUDY)
    threads: 1
    message:
        "DB_generate: {input.sample_metadata} -> {output.sample_info_final}"
    conda:
        "../envs/metagenomics_db.yaml"
    shell:
        """
        if [ -z "{params.study}" ]; then
            python workflow/scripts/metagenomics_db/main.py -d {params.input_dir} -v {params.ver} -o {params.output_dir} -m {input.sample_metadata} -b {params.db_size} &> {log}
        else
            python workflow/scripts/metagenomics_db/main.py -s {params.study} -v {params.ver} -o {params.output_dir} -m {input.sample_metadata} -b {params.db_size} &> {log}
        fi
        """
