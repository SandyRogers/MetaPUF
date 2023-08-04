import os
import pandas as pd

CONTIG_INFO_FILE_DIR = os.path.join(config["outputdir"], "assemblies")

PROCESSED_REPORTS_DIR = os.path.join(config["outputdir"], "Processed_Peptide_Reports")

PROCESSED_RPT = expand("{fname}/processed_{sname}_peptide_report.csv", fname=PROCESSED_REPORTS_DIR, sname=config['samples'])

SAMPLEINFO_FILE    = config["sample_metadata"]
sample_info     = pd.read_csv(SAMPLEINFO_FILE, sep=',')
Assemblies      = sample_info['Assembly'].drop_duplicates().to_list()

GFF_FILE = expand("{fname}/results/{aname}_expressed_proteins.csv", fname=PROCESSED_REPORTS_DIR, aname=Assemblies)

STUDY = config["parameters"]["study"]
PRIDE_ID = config["parameters"]["pride_id"]

SAMPLEINFO_FILE_FINAL = os.path.join(CONTIG_INFO_FILE_DIR, "sample_info_final.csv")

rule gff_format_file:
    input:
        metap_sample_info=SAMPLEINFO_FILE_FINAL,
        rpt=PROCESSED_RPT
    output:
        gff_file=GFF_FILE
    params:
        pride_id=PRIDE_ID,
        reports_dir=PROCESSED_REPORTS_DIR,
        metag_dir=CONTIG_INFO_FILE_DIR
    log:
        expand("logs/{iname}_gff_generate.log", iname=STUDY)
    threads: 1
    conda: "../envs/pandas_mono.yaml"
    message:
        "Generating GFF format file: {input.metap_sample_info} -> {log}"
    shell:
        "python workflow/scripts/gff_generation/main.py -i {input.metap_sample_info} -r {params.reports_dir} -m {params.metag_dir} -p {params.pride_id}"
