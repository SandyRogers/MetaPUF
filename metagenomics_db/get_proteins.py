#-*- coding: utf-8 -*-

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
import gzip
import re
import hashlib
import sys
from os import stat

import pandas as pd
from Bio import SeqIO

def get_fasta(mgnify_study_id: str, fasta_dir: str):
    analyses_request = requests.get(f'https://www.ebi.ac.uk/metagenomics/api/v1/studies/{mgnify_study_id}/analyses')
    analyses = analyses_request.json()
    for analysis in analyses['data']:
        analysis_id = analysis['id']
        downloads_url = analysis['relationships']['downloads']['links']['related']
        downloads_request = requests.get(downloads_url)
        downloads = downloads_request.json()
        for download in downloads['data']:
            if download['attributes']['description']['label'] == 'Predicted CDS (aa)':
                fasta_url = download['links']['self']
                fasta_request = requests.get(fasta_url)
                cds_fasta=os.path.join(fasta_dir,analysis_id+".fasta.gz")
                open(cds_fasta, 'wb').write(fasta_request.content)

def get_protein_metadata(analyses: str, fasta_dir: str, output_dir: str):
    for fasta_in in os.listdir(fasta_dir):
        with gzip.open(os.path.join(output_dir, fasta_in), 'wt') as fasta_out, open(os.path.join(output_dir, analysis_id+"_metadata.txt"), 'w') as peptides_metadata:
            with gzip.open(os.path.join(fasta_dir, fasta_in), 'rt') as cds_fasta:
                for record in SeqIO.parse(cds_fasta, "fasta"):
                    contig_name=record.id
                    sequence = record.seq
                    partial_field, partial, start_coordinate, stop_coordinate, strand, caller = parsing_header(record.description)
                    digest = create_digest(partial + sequence)
                    peptides_metadata.write('\t'.join([analysis_id, contig_name, digest, partial_field, start_coordinate, stop_coordinate, strand, caller]) + '\n')
                    #rename protein in fasta file and save in output_dir
                    record.id=digest
                    SeqIO.write(record, fasta_out, "fasta")


def uniq_proteins(d_dir: str, db_name: str):
    """
    generate the protein search database with unique protein sequences
    :param d_dir:  path of directory containing databases
    :param db_name: name of the search database
    """
    unique_records = {}
    protein_file = os.path.join(d_dir, db_name)
    uniq_db = os.path.join(d_dir, "unique_" + db_name)
    for record in SeqIO.parse(protein_file, "fasta"):
        # if str(record.id) not in unique_records.keys():
        if str(record.id) not in unique_records:
            unique_records[str(record.id)] = record.seq
    with open(uniq_db, "w") as fout:
        for k, v in unique_records.items():
            fout.write(">" + k + "\n")
            fout.write(str(v) + "\n")

def remove_file(file_name):
    try:
        os.remove(file_name)
    except OSError as e:  ## if failed, report it back to the user ##
        print ("Error: %s - %s." % (e.filename, e.strerror))
