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
import time
import urllib.request as request
from contextlib import closing
from os import stat
from traceback import format_exc
from argparse import ArgumentParser

import requests

def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)


def main():  # noqa: C901
    """
    Generates the assembly metadata file with sample info
    """
    parser = ArgumentParser(
        description="Generate protein sequence database for Metagenomics and /or Metatranscriptomics datasets"
    )
    parser.add_argument(
        "-s","--study",
        type=str,
        required=True,
        help="Secondary study accession of the assembled study starting with ERP/SRP/DRP ",
    )
    parser.add_argument(
        "-i","--input_dir",
        type=dir_path,
        required=True,
        help="full path of the study folder containing raw folder",
    )
    starttime = time.time()
    args = parser.parse_args()
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
        + ena_portal_query.format(args.study)
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
    study_file = os.path.join(args.input_dir, args.study + "_info.tsv")
    with open(study_file, "wb") as fout:
        fout.write(file_data)
    logging.info("Metadata file downloaded")
    logging.info("Runtime is {} seconds".format(time.time() - starttime))


if __name__ == "__main__":
    log_file = "assembly_metadata.log"
    logging.basicConfig(
        level=logging.DEBUG, filemode="w", format="%(message)s", datefmt="%H:%M:%S"
    )
    main()
