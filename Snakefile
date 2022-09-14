import os
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
TMPDIR      = os.environ.get("TMPDIR", config['tmp_dir'])
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
# I_PROTEINS = os.environ.get("PROTEINS", config["parameters"]["Proteins"])
SEARCHGUI_OUTPUT = os.environ.get("SEARCHGUI_OUTPUT", config["output"]["searchgui_folder"])
PEPTIDERSHAKER_OUTPUT = os.environ.get("PEPTIDERSHAKER_OUTPUT", config["output"]["peptidershaker_folder"])

##################################################
# WORKDIR
##################################################
workdir:
    WORKDIR



#input files
SAMPLEINFO_FILE = os.path.join(CONFIGDIR, METADATA)
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
PEPTIDESHAKER_MZID = expand("peptideshaker/{fname}_peptideshaker.mzid", fname=Samples)
# FINAL_MZID = expand("{fname}/{fname}_final.mzid", fname=THERMOFOLD)
# RPT_NAMES = [os.path.splitext(os.path.basename(f))[0] for f in Samples]
PSM_RPT = expand("results/reports/proteins/{fname}_psm_report.txt", fname=Samples)
PROTEIN_RPT = expand("results/reports/proteins/{fname}_protein_report.txt", fname=Samples)
PEPTIDE_RPT = expand("results/reports/peptides/{fname}_peptide_report.txt", fname=Samples)
PROCESSED_RPT = expand("PROCESSED_REPORTS_DIR/processed_{fname}_peptide_report.txt", fname=Samples)

# tools
THERMO_EXE = os.path.join(BINDIR, "ThermoRawFileParser/ThermoRawFileParser.exe")
SEARCHGUI_JAR = os.path.join(BINDIR, "SearchGUI-4.0.41/SearchGUI-4.0.41.jar")
SEARCHGUI_PAR_PARAMS = " ".join(["-%s %s" % (k, "'%s'" % v if isinstance(v, str) else str(v)) for k, v in config["searchgui"]["par"].items()])
PEPTIDESHAKER_JAR = os.path.join(BINDIR, "PeptideShaker-2.0.33/PeptideShaker-2.0.33.jar")

#output files
CLUSTER_REPORT = os.path.join(CONTIG_INFO_FILE_DIR, "databases/cluster_report.txt")
PROTEINS_DECOY = os.path.join(TMPDIR, "proteins_decoy_generated_check.txt")
SAMPLEINFO_FILE_FINAL = os.path.join(CONTIG_INFO_FILE_DIR, "sample_info_final.csv")
GFF_FILE_H = expand("PROCESSED_REPORTS_DIR/results/{aname}_high_confidence.gff",aname=Assemblies)
GFF_FILE_L = expand("PROCESSED_REPORTS_DIR/results/{aname}_low_confidence.gff",aname=Assemblies)

##################################################
# RULES
# Each subsequent output file needs to have its target path specified at the beginning.
##################################################
rule ALL:
    input:
        # dynamic(expand("assemblies/databases/unique_{iname}_cluster_set_{{PART}}.faa", iname=STUDY)),
        database=[OUTPUT_FILE, PROTEIN_FILE, SAMPLEINFO_FILE_FINAL, CLUSTER_REPORT]
        # thermo=[THERMORAW, THERMOMGF],
        # searchgui=[PROTEINS_DECOY, SEARCHGUI_PAR, SEARCHGUI_ZIP],
        # report=[PROTEIN_RPT, PEPTIDE_RPT],
        # peptideshaker=PEPTIDESHAKER_MZID,
        # # assembly_list=ASSEMBLY_NAMES,
        # processed=PROCESSED_RPT,
        # gff_files=[GFF_FILE_H, GFF_FILE_L]



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
            shell("python metagenomics_db/main.py -s {params.study} -v {params.ver} -o {params.output_dir} -m {input.sample_metadata} -b {params.db_size} ")
        elif params.study in config["parameters"]["Input_dir"]:
            shell("python metagenomics_db/main.py -i {params.study} -v {params.ver} -o {params.output_dir} -m {input.sample_metadata} -b {params.db_size} ")


#########################
# ThermoRawFileParser
#########################
# https://github.com/compomics/ThermoRawFileParser
# rule fetch_raw_files:
#     input:
#         RawURLs
#     output:
#         THERMORAW
#     log:
#         expand("logs/{fname}_fetch_raw_files.log",fname=PRIDE_ID)
#     shell:
#         "wget -P input/Raw -i {input}"


rule thermorawfileparser:
    input:
        exe=THERMO_EXE,
        info=SAMPLEINFO_FILE
    output:
        raw=THERMORAW,
        mgf=THERMOMGF
    params:
        folder=THERMOFOLD
    # log:
        # expand("logs/{fname}_thermorawfileparser.log",fname=PRIDE_ID)
    threads: 1
    message:
        "ThermoRawFileParser: {input.info} -> {output.mgf}"
    shell:
        "python coping_raw_files.py -exe {input.exe} -info {input.info} -out {params.folder}"
    # shell:
    #     "mono {input.exe} -d=$(dirname {input.raws[0]}) -o=$(dirname {output[0]}) -f=1 -m=0 &> {log}"


#########################
# SearchGUI
#########################
rule searchgui_decoy:
    input:
        info=SAMPLEINFO_FILE_FINAL,
        human_db=HUMAN_FASTA,
        crap_db=CRAP_FASTA,
        jar=SEARCHGUI_JAR,
        cluster_rpt=CLUSTER_REPORT
    output:
        PROTEINS_DECOY
    log:
        expand("logs/{fname}_SearchGUI_decoy.log",fname=PRIDE_ID)
    params:
        tmpdir = TMPDIR
        # logdir = "logs/SearchGUI_decoy"
    threads: 1
    conda:
        os.path.join(ENVDIR, "IMP_proteomics.yaml")
    message:
        "SearchGUI decoy: {input.cluster_rpt} -> {output}"
    shell:
        "python generating_decoy.py -jar {input.jar} -info {input.info} "
        "-human {input.human_db} -crap {input.crap_db} -tmp {params.tmpdir}"
    # shell:
    #     "for protein in {input.faa}; do cat {input.human_db} {input.crap_db} >> $protein; "
    #     "java -cp {input.jar} eu.isas.searchgui.cmd.FastaCLI -in $protein "
    #     "-decoy -temp_folder {params.tmpdir} -log {params.logdir} &> {log}; done "


rule searchgui_config:
    input:
        jar=SEARCHGUI_JAR
    output:
        SEARCHGUI_PAR
    log:
        expand("logs/{fname}_SearchGUI_params.log",fname=PRIDE_ID)
    params:
        params = SEARCHGUI_PAR_PARAMS,
        tmpdir = TMPDIR
        # logdir = "logs/SearchGUI_params"
    threads: 1
    conda:
        os.path.join(ENVDIR, "IMP_proteomics.yaml")
    message:
        "SearchGUI parameters: {input} -> {output}"
    shell:
        "java -cp {input.jar} eu.isas.searchgui.cmd.IdentificationParametersCLI -out {output} "
        "{params.params} -temp_folder {params.tmpdir} &> {log}"


rule searchgui_search:
    input:
        par=SEARCHGUI_PAR,
        flag=PROTEINS_DECOY,
        # mgf=THERMOMGF,
        jar=SEARCHGUI_JAR,
        info=SAMPLEINFO_FILE
    output:
        SEARCHGUI_ZIP
    # log:
        # expand("logs/{fname}_SearchGUI_search.log",fname=PRIDE_ID)
    params:
        # name=expand("{fname}_searchgui", fname=PRIDE_ID),
        # tmpdir = TMPDIR,
        logdir = "logs/SearchGUI_search"
    threads: 10
    conda:
        os.path.join(ENVDIR, "IMP_proteomics.yaml")
    message:
        "SearchGUI search: {input.flag} -> {output}"
    shell:
        "python searchgui_search.py -s -jar {input.jar} -in {input.info} "
        "-out $(dirname {output[0]}) -par {input.par}"
    # shell:
    #     """
    #     java -cp {input.jar} eu.isas.searchgui.cmd.SearchCLI \
    #         -spectrum_files $(dirname {input.mgf[0]}) \
    #         -fasta_file {input.faa} \
    #         -output_folder $(dirname {output}) \
    #         -id_params {input.par} \
    #         -xtandem 1 \
    #         -msgf 1 \
    #         -comet 0 \
    #         -andromeda 0 \
    #         -threads {threads} \
    #         -output_default_name {params.name} \
    #         -output_option 0 \
    #         -output_data 1 \
    #         -output_date 0 \
    #         -log {params.logdir} \
    #         &> {log} && touch {output}
    #     """


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
        outputdir=PEPTIDERSHAKER_OUTPUT,
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
    # shell:
    #     "for search in {input.searchgui}; do java -cp {input.jar} eu.isas.peptideshaker.cmd.PeptideShakerCLI "
    #     "-reference 'peptideshaker_peptideshaker_1' "
    #     "-identification_files $search "
    #     "-out_reports $(dirname {output.protein[0]}) -reports 6,9 "
    #     "-report_prefix $(cut -d'_' -f1 <<<${search##*/}) "
    #     "-output_file {output.mzid[0]} -contact_first_name 'Shengbo' -contact_last_name 'Wang' "
    #     "-contact_email 'shengbo_wang@ebi.ac.uk' -contact_address 'EBI' -organization_name 'EBI' "
    #     "-organization_email 'test@ebi.ac.uk' -organization_address 'Cambridge'; "
    #     "-threads {threads} &> {log}; done"



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
        "Post-processing: {input} -> {output}"
    shell:
        "python post_report_generation/main.py -s {input} -p {params} &> {log}"

#########################
# Gff format file
#########################
if  config["parameters"]["Study"]!='':
    rule gff_format_file:
        input:
            metap_sample_info=SAMPLEINFO_FILE_FINAL,
            rpt=PROCESSED_RPT
        output:
            gff_file_high_confidence=GFF_FILE_H,
            gff_file_low_confidence=GFF_FILE_L
        params:
            reports_dir=PROCESSED_REPORTS_DIR,
            metag_dir=CONTIG_INFO_FILE_DIR,
            pride_id=PRIDE_ID
        threads: 1
        message:
            "Generating GFF format file: {input.metap_sample_info} -> {output.gff_file_high_confidence}"
        shell:
            "python gff_generation/main.py -s {input.metap_sample_info} -r {params.reports_dir} "
            "-m {params.metag_dir} -p {params.pride_id} "
else:
    sys.exit()
