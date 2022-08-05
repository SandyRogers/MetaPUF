# -*- coding: utf-8 -*-

# Copyright 2022 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import logging
import os
import shutil
import re
import time
import urllib.request as request
from contextlib import closing
from os import stat
from traceback import format_exc

import requests


def check_study_accession(study_accession):
    study_accession_re = re.compile(r"((E|D|S)RP[0-9]{6,})")
    if not study_accession or not study_accession_re.match(study_accession):
        logging.error(f"{study_accession} is not an valid Secondary study accession.")
        sys.exit(1)
    else:
        print("Starting the workflow")

def assembly_info(study: str, fpath: str):
    """
    Get information about the study of interest from ENA
    :param study: secondary study accession starting ERP
    :param fpath: file path

    :return: a tsv file with all the information from the api
    """
    ena_portal_base_api_url = "https://www.ebi.ac.uk/ena/portal/api/search?"
    ena_portal_params = [
        "dccDataOnly=false",
        "result=analysis",
        "download=True",
        "format=tsv",
    ]
    ena_portal_query = "query=secondary_study_accession=%22{}%22&fields="
    ena_portal_fields = [
        "secondary_study_accession",
        "secondary_sample_accession",
        "sample_alias",
        "analysis_accession",
        "assembly_type",
        "generated_ftp",
    ]
    ena_portal_api_url = (
        ena_portal_base_api_url
        + "&".join(ena_portal_params)
        + "&"
        + ena_portal_query.format(study)
        + ",".join(ena_portal_fields)
    )
    try:
        download = requests.get(ena_portal_api_url)
        file_data = download.content
        download.raise_for_status()
    except requests.exceptions.Timeout:
        time.sleep(10)
    except requests.exceptions.ConnectionError as exception:
        logging.info(f"{ena_portal_base_api_url} : {exception!r} \n {format_exc()} \n")
        raise SystemExit(exception)
    except requests.exceptions.HTTPError as exception:
        logging.info(f"{ena_portal_base_api_url} : {exception!r} \n {format_exc()} \n")
        raise SystemExit(exception)
    study_file = os.path.join(fpath, study + "_info.tsv")
    with open(study_file, "wb") as fout:
        fout.write(file_data)


def ftp_download(fpath: str, ids: str):
    """
    Download assembly files (.fasta.gz) for a study to specified location
    :param fpath: path of project_id legacy_file_location
    :param ids: a text file contating the ftp path for all the .fasta.gz files
    """
    with open(ids) as fin:
        headers = (fin.readline()).strip().split("\t")
        url_index = headers.index("generated_ftp")
        analysis_index = headers.index("analysis_accession")
        for line in fin:
            str1 = line.strip().split("\t")
            url = str1[url_index]
            file_extension = ".".join((url.strip().split("/")[-1]).split(".")[1:])
            base_name = str1[analysis_index] + "." + file_extension
            file_name = os.path.join(fpath, base_name)
            if (not os.path.exists(file_name)) or (stat(file_name).st_size == 0):
                with closing(request.urlopen("ftp://" + url)) as r:
                    logging.info("Downloading file from FTP server..." + url)
                    logging.info("Getting the file...")
                    with open(file_name, "wb") as fout:
                        shutil.copyfileobj(r, fout)
                    logging.info("File " + file_name + " downloaded.")
