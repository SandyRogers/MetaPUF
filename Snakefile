import os
import sys
import pandas as pd

##################################################
# SHELL
##################################################
# default executable for snakmake
shell.executable("bash")

##################################################
# PATHS
##################################################
# default configuration file
configfile:
    srcdir("config/config.proteomics.yaml")

# relevant paths
BINDIR      = srcdir("workflow/bin")
ENVDIR      = srcdir("workflow/envs")
CONFIGDIR   = srcdir("config")
RESOURCEDIR = srcdir("resources")
WORKDIR     = os.environ.get("WORKDIR", config['workdir'])
OUTPUTDIR   = os.environ.get("OUTPUTDIR", config['outputdir'])
TMPDIR      = os.environ.get("TMPDIR", config['tmpdir'])
CONTIG_INFO_FILE_DIR = os.path.join(OUTPUTDIR,"assemblies")
PROCESSED_REPORTS_DIR = os.path.join(OUTPUTDIR,"Processed_Peptide_Reports")


# parameters
FIRST_NAME  = os.environ.get("FIRST_NAME", config["parameters"]["first_name"])
LAST_NAME   = os.environ.get("LAST_NAME", config["parameters"]["last_name"])
EMAIL       = os.environ.get("EMAIL", config["parameters"]["email"])
ADDRESS     = os.environ.get("ADDRESS", config["parameters"]["address"])
ORG_NAME    = os.environ.get("ORG_NAME", config["parameters"]["org_name"])
ORG_EMAIL   = os.environ.get("ORG_EMAIL", config["parameters"]["org_email"])
ORG_ADDRESS = os.environ.get("ORG_ADDRESS", config["parameters"]["org_address"])
THERMOFOLD  = os.environ.get("THERMOFOLD", config["parameters"]["ThermoFold"])
STUDY       = os.environ.get("STUDY", config["parameters"]["Study"])
PRIDE_ID    = os.environ.get("PRIDE_ID", config["parameters"]["Pride_id"])
VERSION     = os.environ.get("VERSION", config["parameters"]["Version"])
DB_SIZE     = os.environ.get("DB_SIZE", config["parameters"]["Db_size"])
METADATA    = os.environ.get("METADATA", config["raws"]["Metadata"])
SEARCHGUI_OUTPUT = os.environ.get("SEARCHGUI_OUTPUT", config["output"]["searchgui_folder"])
PEPTIDESHAKER_OUTPUT = os.environ.get("PEPTIDESHAKER_OUTPUT", config["output"]["peptideshaker_folder"])


##################################################
# WORKDIR
##################################################
workdir:
    OUTPUTDIR

#input files
SAMPLEINFO_FILE = os.path.join(CONFIGDIR, METADATA)
CONFIG_FILE     = os.path.join(CONFIGDIR, "config.proteomics.yaml")
sample_info     = pd.read_csv(SAMPLEINFO_FILE, sep=',')
Assemblies      = sample_info['Assembly'].drop_duplicates().to_list()
Samples         = sample_info['Sample'].drop_duplicates().to_list()
sample_raw      = sample_info[['Sample','Raw file']].drop_duplicates()
sample_raw[['filename','extension']] = sample_raw['Raw file'].str.split('.',expand=True)
RawFileNames    = (sample_raw['Sample']+'/'+sample_raw['filename']).to_list()
CRAP_FASTA      = os.path.join(RESOURCEDIR, "crap_db.fa")
HUMAN_FASTA     = os.path.join(RESOURCEDIR, "human_db.fa")

# data (output from one rule being used as input for another))
OUTPUT_FILE = expand("assemblies/{aname}_contig_info.txt", aname=Assemblies)
PROTEIN_FILE = expand("assemblies/{aname}.faa.gz", aname=Assemblies)
THERMORAW = expand("{iname}/{bname}.raw", iname=THERMOFOLD, bname=RawFileNames)
THERMOMGF = expand("{iname}/{bname}.mgf", iname=THERMOFOLD, bname=RawFileNames)
SEARCHGUI_PAR  = expand("searchgui/{fname}_searchgui.par", fname=PRIDE_ID)
SEARCHGUI_ZIP  = expand("searchgui/{fname}_searchgui.zip", fname=Samples)
PEPTIDESHAKER_MZID = expand("{fname}/{sname}_peptideshaker.mzid", fname=PEPTIDESHAKER_OUTPUT, sname=Samples)
PSM_RPT = expand("results/reports/proteins/{fname}_psm_report.txt", fname=Samples)
PROTEIN_RPT = expand("{fname}/peptideshaker_{sname}_1_Default_Protein_Report.txt", fname=PEPTIDESHAKER_OUTPUT, sname=Samples)
PEPTIDE_RPT = expand("{fname}/peptideshaker_{sname}_1_Default_Peptide_Report.txt", fname=PEPTIDESHAKER_OUTPUT, sname=Samples)
PROCESSED_RPT = expand("{fname}/processed_{sname}_peptide_report.csv", fname=PROCESSED_REPORTS_DIR, sname=Samples)

# tools
THERMO_EXE = os.path.join(BINDIR, "ThermoRawFileParser/ThermoRawFileParser.exe")
SEARCHGUI_JAR = os.path.join(BINDIR, "SearchGUI-3.3.20/SearchGUI-3.3.20.jar")
SEARCHGUI_PAR_PARAMS = " ".join(["-%s %s" % (k, "'%s'" % v if isinstance(v, str) else str(v)) for k, v in config["searchgui"]["par"].items()])
PEPTIDESHAKER_JAR = os.path.join(BINDIR, "PeptideShaker-1.15.46/PeptideShaker-1.15.46.jar")

#output files
CLUSTER_REPORT = os.path.join(CONTIG_INFO_FILE_DIR, "databases/cluster_report.txt")
PROTEINS_CHECK = os.path.join(CONTIG_INFO_FILE_DIR, "databases/proteins_decoy_params_generated_check.txt")
SAMPLEINFO_FILE_FINAL = os.path.join(CONTIG_INFO_FILE_DIR, "sample_info_final.csv")
GFF_FILE = expand("{fname}/results/{aname}_expressed_proteins.csv", fname=PROCESSED_REPORTS_DIR, aname=Assemblies)

##################################################
# RULES
# Each subsequent output file needs to have its target path specified at the beginning.
##################################################
rule ALL:
    input:
        database=[OUTPUT_FILE, PROTEIN_FILE, SAMPLEINFO_FILE_FINAL, CLUSTER_REPORT],
        thermo=THERMOMGF,
        searchgui=[PROTEINS_CHECK, SEARCHGUI_ZIP],
        report=[PROTEIN_RPT, PEPTIDE_RPT],
        peptideshaker=PEPTIDESHAKER_MZID,
        processed=PROCESSED_RPT,
        gff_files=GFF_FILE


#########################
# Generate protein search database
#########################
rule generate_db:
    input:
        sample_metadata=SAMPLEINFO_FILE
    output:
        contigs_dir=OUTPUT_FILE,
        protein_file=PROTEIN_FILE,
        sample_info_final=SAMPLEINFO_FILE_FINAL,
        cluster_rpt=CLUSTER_REPORT
    params:
        study=STUDY if config["parameters"]["Study"] else config["parameters"]["Input_dir"],
        ver=VERSION,
        output_dir=OUTPUTDIR,
        db_size=DB_SIZE
    log:
        expand("logs/{iname}_db_generate.log",iname=STUDY)
    threads: 1
    message:
        "DB_generate: {input.sample_metadata} -> {output.sample_info_final}"
    run:
        if params.study in config["parameters"]["Study"]:
            shell("python metagenomics_db/main.py -s {params.study} -v {params.ver} -o {params.output_dir} -m {input.sample_metadata} -b {params.db_size} &> {log}")
        elif params.study in config["parameters"]["Input_dir"]:
            shell("python metagenomics_db/main.py -i {params.study} -v {params.ver} -o {params.output_dir} -m {input.sample_metadata} -b {params.db_size} &> {log}")


#########################
# ThermoRawFileParser
#########################
# https://github.com/compomics/ThermoRawFileParser
rule thermorawfileparser:
    input:
        exe=THERMO_EXE,
        info=SAMPLEINFO_FILE
    output:
        mgf=THERMOMGF
    params:
        folder=THERMOFOLD
    threads: 1
    message:
        "ThermoRawFileParser: {input.info} -> {output.mgf}"
    shell:
        "python coping_raw_files.py -exe {input.exe} -info {input.info} -out {params.folder}"


#########################
# SearchGUI
#########################
rule searchgui_decoy:
    input:
        conf=CONFIG_FILE,
        info=SAMPLEINFO_FILE_FINAL,
        human_db=HUMAN_FASTA,
        crap_db=CRAP_FASTA,
        jar=SEARCHGUI_JAR,
        cluster_rpt=CLUSTER_REPORT
    output:
        PROTEINS_CHECK
    threads: 1
    conda:
        os.path.join(ENVDIR, "IMP_proteomics.yaml")
    message:
        "SearchGUI decoy: {input.cluster_rpt} -> {output}"
    shell:
        "python generating_decoy.py -jar {input.jar} -info {input.info} "
        "-human {input.human_db} -crap {input.crap_db} -p {input.conf}"


rule searchgui_search:
    input:
        flag=PROTEINS_CHECK,
        mgf=THERMOMGF,
        jar=SEARCHGUI_JAR,
        info=SAMPLEINFO_FILE_FINAL
    output:
        SEARCHGUI_ZIP
    threads: 10
    conda:
        os.path.join(ENVDIR, "IMP_proteomics.yaml")
    message:
        "SearchGUI search: {input.mgf} -> {output}"
    shell:
        "python searchgui_search.py -s -jar {input.jar} -in {input.info} "
        "-out $(dirname {output[0]}) "


#########################
# PeptideShaker
#########################
# http://compomics.github.io/projects/peptide-shaker
rule peptideshaker_load:
    input:
        searchgui=SEARCHGUI_ZIP,
        jar=PEPTIDESHAKER_JAR,
        info=SAMPLEINFO_FILE
    output:
        protein=PROTEIN_RPT,
        peptide=PEPTIDE_RPT,
        mzid=PEPTIDESHAKER_MZID
    params:
        outputdir=PEPTIDESHAKER_OUTPUT,
        fn=FIRST_NAME,
        ln=LAST_NAME,
        ce=EMAIL,
        ca=ADDRESS,
        on=ORG_NAME,
        oe=ORG_EMAIL,
        oa=ORG_ADDRESS
    log:
        expand("logs/{fname}_PeptideShaker_load.log",fname=Samples)
    threads: 10
    conda:
        os.path.join(ENVDIR, "IMP_proteomics.yaml")
    message:
        "PeptideShaker load SearchGUI results: {input.info} -> {output.mzid}, {output.protein}, {output.peptide}"
    shell:
        "python searchgui_search.py -p -jar {input.jar} -in {input.info} "
        "-out {params.outputdir} -fn {params.fn} -ln {params.ln} -ce {params.ce} -ca {params.ca} "
        "-on {params.on} -oe {params.oe} -oa {params.oa}"


########################
# Generate post processing reports
########################
rule post_processing:
    input:
        info=SAMPLEINFO_FILE,
        protein=PROTEIN_RPT,
        peptide=PEPTIDE_RPT
    output:
        PROCESSED_RPT
    params:
        PRIDE_ID
    log:
        expand("logs/{fname}_post_processing.log", fname=PRIDE_ID)
    threads: 1
    message:
        "Post-processing: {input.info} -> {output}"
    shell:
        "python post_report_generation/main.py -s {input.info} -p {params} &> {log}"


#########################
# Gff format file
#########################
rule gff_format_file:
    input:
        metap_sample_info=SAMPLEINFO_FILE_FINAL,
        rpt=PROCESSED_RPT
    output:
        gff_file=GFF_FILE
    params:
        study=STUDY if config["parameters"]["Study"] else config["parameters"]["Input_dir"],
        pride_id=PRIDE_ID,
        reports_dir=PROCESSED_REPORTS_DIR,
        metag_dir=CONTIG_INFO_FILE_DIR
    log:
        expand("logs/{iname}_gff_generate.log",iname=STUDY)
    threads: 1
    message:
        "Generating GFF format file: {input.metap_sample_info} -> {log}"
    run:
        if params.study in config["parameters"]["Study"]:
            shell("python gff_generation/main.py -s {params.study} -i {input.metap_sample_info} -r {params.reports_dir} -m {params.metag_dir} -p {params.pride_id} ")
        elif params.study in config["parameters"]["Input_dir"]:
            shell("python gff_generation/main.py -d {params.study} -i {input.metap_sample_info} -r {params.reports_dir} -m {params.metag_dir} -p {params.pride_id} ")
