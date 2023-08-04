import os

SAMPLEINFO_FILE    = config["sample_metadata"]
THERMOFOLD = os.path.join(config["outputdir"], "thermofold")

THERMOMGF = expand("{iname}/{bname}.mgf", iname=THERMOFOLD, bname=config['raw_file_names'])

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

thermorawfileparser_path = expand(
    "{tooldir}/{exe}",
    tooldir = config['thermo']['tooldir'],
    exe = config['thermo']['exe']
)

thermorawfileparser_tool = config['thermo']['tooldir']

rule download_thermorawfileparser:
    input:
        exe_zip=HTTP.remote(expand(config['thermo']['url'], version=config['thermo']['version']), keep_local=False)
    output:
        exe=thermorawfileparser_path
    params:
        thermorawfileparser_tool=thermorawfileparser_tool
    log: "logs/download_thermorawfileparser.log"
    shell:
        "unzip {input.exe_zip} -d {params.thermorawfileparser_tool} > {log}"

rule thermorawfileparser:
    input:
        info=SAMPLEINFO_FILE,
        exe=thermorawfileparser_path
    output:
        mgf=THERMOMGF,
        folder=directory(THERMOFOLD)
    params:
        folder=THERMOFOLD,
        raw=config["parameters"]["raw_dir"]
    threads: 1
    conda:
        "../envs/pandas_mono.yaml"
    message:
        "ThermoRawFileParser: {input.info} -> {output.mgf}"
    shell:
        "python workflow/scripts/coping_raw_files.py -exe {input.exe} -info {input.info} -in {params.raw} -out {params.folder}"
