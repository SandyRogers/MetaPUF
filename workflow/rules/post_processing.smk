import os

SAMPLEINFO_FILE    = config["sample_metadata"]
PROCESSED_REPORTS_DIR = os.path.join(config["outputdir"], "Processed_Peptide_Reports")

PEPTIDESHAKER_OUTPUT = os.path.join(config["outputdir"], "peptideshaker")
PROTEIN_RPT = expand("{fname}/peptideshaker_{sname}_1_Default_Protein_Report.txt", fname=PEPTIDESHAKER_OUTPUT, sname=config['samples'])
PEPTIDE_RPT = expand("{fname}/peptideshaker_{sname}_1_Default_Peptide_Report.txt", fname=PEPTIDESHAKER_OUTPUT, sname=config['samples'])
PROCESSED_RPT = expand("{fname}/processed_{sname}_peptide_report.csv", fname=PROCESSED_REPORTS_DIR, sname=config['samples'])
PRIDE_ID = config["parameters"]["pride_id"]

rule post_processing:
    input:
        info=SAMPLEINFO_FILE,
        protein=PROTEIN_RPT,
        peptide=PEPTIDE_RPT
    output:
        PROCESSED_RPT
    params:
        pid=PRIDE_ID,
        dir=PROCESSED_REPORTS_DIR
    log:
        expand("logs/{fname}_post_processing.log", fname=PRIDE_ID)
    threads: 1
    conda: "../envs/pandas_mono.yaml"
    message:
        "Post-processing: {input.info} -> {output}"
    shell:
        "python workflow/scripts/post_report_generation/main.py -s {input.info} -d {params.dir} -p {params.pid} &> {log}"
