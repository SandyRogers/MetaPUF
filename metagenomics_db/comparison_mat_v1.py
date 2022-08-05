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

import sys
import logging
import re
import subprocess


def check_study_accession(study_accession):
    study_accession_re = re.compile(r"((E|D|S)RP[0-9]{6,})")
    if not study_accession or not study_accession_re.match(study_accession):
        logging.error(f"{study_accession} is not an valid Secondary study accession.")
        sys.exit(1)
    else:
        print("Starting the workflow")

def sourmash_sig(in_file: str, out_file:str):
    """
    Compute sourmash signatures for the .fasta files
    :param in_file: path of the  input file
    :param out_file: path of the  output file

    """
    if not in_file or not isinstance(in_file, str):
        raise ValueError("Fastafiles is not valid or empty!")
    cmd_sig = "  ".join(["sourmash sketch dna -p  k=31,scaled=1000",  in_file, "-o ", out_file])
    print(cmd_sig)
    subprocess.call(cmd_sig, shell=True)

def signature_compare(file_out: str, matrix_file: str, sig_dir: str):
    """
    Compute similarity matrix
    :param sig_dir: path of the signature directories
    :param file_out: name of the study
    :param matrix_file: name of the file for saving the comparidon matrix
    """
    if not sig_dir or not isinstance(sig_dir, str):
        raise ValueError("Sourmash signatures directory is empty!")
    cmd_index = "  ".join(["sourmash compare -k 31",sig_dir+"/*.sig   --csv", file_out, "-o ", matrix_file])
    subprocess.call(cmd_index, shell=True)
