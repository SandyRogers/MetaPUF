import os

SAMPLEINFO_FILE    = config["sample_metadata"]
PEPTIDESHAKER_OUTPUT = os.path.join(config["outputdir"], "peptideshaker")
PEPTIDESHAKER_MZID = expand("{fname}/{sname}_peptideshaker.mzid", fname=PEPTIDESHAKER_OUTPUT, sname=config['samples'])
PROTEIN_RPT = expand("{fname}/peptideshaker_{sname}_1_Default_Protein_Report.txt", fname=PEPTIDESHAKER_OUTPUT, sname=config['samples'])
PEPTIDE_RPT = expand("{fname}/peptideshaker_{sname}_1_Default_Peptide_Report.txt", fname=PEPTIDESHAKER_OUTPUT, sname=config['samples'])
SEARCHGUI_ZIP  = expand("{output}/searchgui/{fname}_searchgui.zip", output=config["outputdir"], fname=config['samples'])

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

peptideshaker_jarname = expand(
    config["peptideshaker"]["jar"],
    version = config["peptideshaker"]["version"]
)

peptideshaker_tooldir = config['peptideshaker']['tooldir']

peptideshaker_path = expand(
    "{tooldir}/{jar}",
    tooldir = peptideshaker_tooldir,
    jar = peptideshaker_jarname
)

rule download_peptideshaker:
    input:
        zip=HTTP.remote(expand(config['peptideshaker']['url'], version=config['peptideshaker']['version']), keep_local=False)
    output:
        jar=peptideshaker_path
    params:
        peptideshaker_tooldir=peptideshaker_tooldir
    log: "logs/download_peptideshaker.log"
    shell:
        "unzip {input.zip} -d {params.peptideshaker_tooldir} > {log}"

rule peptideshaker_load:
    input:
        searchgui=SEARCHGUI_ZIP,
        jar=peptideshaker_path,
        info=SAMPLEINFO_FILE
    output:
        protein=PROTEIN_RPT,
        peptide=PEPTIDE_RPT,
        mzid=PEPTIDESHAKER_MZID
    params:
        outputdir=config["outputdir"],
    log:
        expand("logs/{fname}_PeptideShaker_load.log", fname=config['samples'])
    threads: 10
    conda:
        "../envs/pandas_java.yaml"
    message:
        "PeptideShaker load SearchGUI results: {input.info} -> {output.mzid}, {output.protein}, {output.peptide}"
    shell:
        "python workflow/scripts/searchgui_search.py -p -jar {input.jar} -in {input.info} "
        "-out {params.outputdir}"
