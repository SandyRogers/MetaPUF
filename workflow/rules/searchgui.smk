import json
import os

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

CONTIG_INFO_FILE_DIR = os.path.join(config["outputdir"], "assemblies")
CLUSTER_REPORT = os.path.join(CONTIG_INFO_FILE_DIR, "databases/cluster_report.txt")
PROTEINS_CHECK = os.path.join(CONTIG_INFO_FILE_DIR, "databases/proteins_decoy_params_generated_check.txt")
SAMPLEINFO_FILE_FINAL = os.path.join(CONTIG_INFO_FILE_DIR, "sample_info_final.csv")

SEARCHGUI_ZIP  = expand("{output}/searchgui/{fname}_searchgui.zip", output=config["outputdir"], fname=config['samples'])

searchgui_jarname = expand(
    config["searchgui"]["jar"],
    version = config["searchgui"]["version"]
)

searchgui_tooldir = config['searchgui']['tooldir']

searchgui_path = expand(
    "{tooldir}/{jar}",
    tooldir = searchgui_tooldir,
    jar = searchgui_jarname
)

rule download_searchgui:
    input:
        zip=HTTP.remote(expand(config['searchgui']['url'], version=config['searchgui']['version']), keep_local=False)
    output:
        jar=searchgui_path
    params:
        searchgui_tooldir=searchgui_tooldir
    log: "logs/download_searchgui.log"
    shell:
        "tar -xzvf {input.zip} -C {params.searchgui_tooldir} > {log}"


rule searchgui_decoy:
    input:
        info=SAMPLEINFO_FILE_FINAL,
        human_db=config["searchgui"]["resources"]["human_db"],
        crap_db=config["searchgui"]["resources"]["crap_db"],
        jar=searchgui_path,
        cluster_rpt=CLUSTER_REPORT,
    output:
        PROTEINS_CHECK
    params:
        searchgui_params = json.dumps(config["searchgui"]["params"]),
        assemblies_dir = os.path.join(config["outputdir"], "assemblies")
    threads: 1
    conda:
        "../envs/pandas_java.yaml"
    message:
        "SearchGUI decoy: {input.cluster_rpt} -> {output}"
    shell:
        "python workflow/scripts/generating_decoy.py -jar {input.jar} -info {input.info} "
        "-human {input.human_db} -crap {input.crap_db} -sgpars '{params.searchgui_params}' "
        "-adir {params.assemblies_dir}"


THERMOFOLD = os.path.join(config["outputdir"], "thermofold")
THERMOMGF = expand("{iname}/{bname}.mgf", iname=THERMOFOLD, bname=config['raw_file_names'])


rule searchgui_search:
    input:
        flag=PROTEINS_CHECK,
        mgf=THERMOMGF,
        jar=searchgui_path,
        info=SAMPLEINFO_FILE_FINAL,
        raw=THERMOFOLD
    output:
        SEARCHGUI_ZIP
    threads: 10
    params:
        outdir=config["outputdir"]
    conda:
        "../envs/pandas_java.yaml"
    message:
        "SearchGUI search: {input.mgf} -> {output}"
    shell:
        "python workflow/scripts/searchgui_search.py -s -jar {input.jar} -in {input.info} "
        "-raw {input.raw} -out {params.outdir}"
