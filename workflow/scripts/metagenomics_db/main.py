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

import gzip
import logging
import os
import shutil
import subprocess
import sys
import time
from argparse import ArgumentParser
from collections import defaultdict

import fetch_data as fd
import get_clusters as gc
import pandas as pd
import screed
import sourmash


def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)


def main():  # noqa: C901
    """
    Generates protein databases containing non-redundant proteins from
    metagenomics and/or metatranscriptomics assemblies
    """
    parser = ArgumentParser(
        description="Generate protein sequence database for Metagenomics and /or Metatranscriptomics datasets"
    )
    """
    The user can enter their own data using --input_dir option
    or
    select --study to analyse study published on MGnify website
    """

    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-s",
        "--study",
        type=str,
        help="Secondary study accession of the assembled study in MGnify starting with ERP/SRP/DRP ",
    )
    group.add_argument(
        "-d",
        "--input_dir",
        type=dir_path,
        help="full path of the study folder containing protein files ",
    )
    parser.add_argument(
        "-v",
        "--ver",
        type=str,
        required=True,
        help="pipeline version of MGnify analysis for the study",
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        type=dir_path,
        required=True,
        help="full path of the study folder containing results",
    )
    parser.add_argument(
        "-m",
        "--metadata",
        type=str,
        required=True,
        help="full path of the sample-assembly mapping file (.csv) with filename.",
    )
    parser.add_argument(
        "-b",
        "--db_size",
        type=int,
        help="input the maximum size of the protein search database in bytes",
        default=1073741824,
    )

    starttime = time.time()
    args = parser.parse_args()
    sample_assembly_map = defaultdict(set)

    assembly_folder = os.path.join(args.output_dir, "assemblies")
    os.makedirs(assembly_folder, exist_ok=True)

    logging.info(f"Will read samples info from {args.metadata}")
    samples = pd.read_csv(args.metadata, sep=",")

    os.chdir(args.output_dir)
    if args.study:
        logging.info(f"Will download data from MGnify API for study {args.study}")
        fd.check_study_accession(args.study)
        # getting the sequnece data and predicted cds
        cmd_get_data = "  ".join(
            [
                "mg-toolkit -d bulk_download -a",
                args.study,
                "-p ",
                args.ver,
                "-g sequence_data",
            ]
        )
        subprocess.call(cmd_get_data, shell=True)
        sequence_dir = (
            args.output_dir + "/" + args.study + "/" + args.ver + "/sequence_data"
        )

    elif args.input_dir:
        logging.info(f"Will use MGnify data from input dir {args.input_dir}")
        ## Assembly input files should be gzipped and end with "_FASTA.fasta.gz"
        ## Prodigal (v.2.6.3) predicted CDS file should be gzipped and end with "_FASTA_predicted_cds.faa.gz"
        if len(os.listdir(args.input_dir)) == 0:
            sys.exit("{} is empty".format(args.input_dir))
        else:
            sequence_dir = str(args.input_dir)
        for input_file in os.listdir(args.input_dir):
            if input_file.endswith("_FASTA.fasta.gz"):
                logging.info(f"Input file {input_file} is correctly named")
                continue
            elif input_file.endswith("_FASTA_predicted_cds.faa.gz"):
                logging.info(f"Input file {input_file} is correctly named")
                continue
            else:
                logging.warning(
                    f"The input file {input_file} is not named with correct naming convention. Please check documentation for correct naming convention"
                )

    for idx, row in samples.iterrows():
        sample_assembly_map[row["Secondary Sample Accession"]].add(row["Assembly"])
    logging.info("Mapping between samples and the assemblies: ", sample_assembly_map)
    for k, v in sample_assembly_map.items():
        sample_file = os.path.join(assembly_folder, k + ".fasta.gz")
        logging.info(f"Reading assemblies from sample {sample_file}")
        with gzip.open(sample_file, "wt") as wf:
            for assembly in v:
                # renaming headers of the predicted cds in assemblies
                fd.rename_contigs(assembly, assembly_folder, sequence_dir)
                with gzip.open(
                    os.path.join(sequence_dir, assembly + "_FASTA.fasta.gz"), "rt"
                ) as infile:
                    shutil.copyfileobj(infile, wf)
    matrix_file = os.path.join(assembly_folder, "matrix_file.tsv")

    meta_genomes = [
        assembly_folder + "/" + k + ".fasta.gz" for k in sample_assembly_map.keys()
    ]
    minhashes = []
    for contigs in meta_genomes:
        mh = sourmash.MinHash(ksize=31, n=0, scaled=1000)
        for record in screed.open(contigs):
            query_seq = record.sequence
            mh.add_sequence(query_seq, True)
        minhashes.append(mh)
    with open(matrix_file, "w") as f_out:
        for i, e in enumerate(minhashes):
            basename = os.path.basename(meta_genomes[i])
            name = ".".join((basename).split(".")[:-2])
            f_out.write(name + ",")
            for j, e2 in enumerate(minhashes):
                x = e.jaccard(minhashes[j])
                f_out.write(str(round(x, 3)) + ",")
            f_out.write("\n")
    data = pd.read_csv(matrix_file, sep=",", header=None, index_col=[0])
    data = data.dropna(axis="columns", how="all")
    data.index.names = ["index"]
    col_names = data.index.to_list()
    data.columns = col_names
    logging.info(f"matrix data frame: {data}")
    # returns a dictionary containing assembly name and count of proteins in the assembly

    proteins_info = fd.count_proteins(sample_assembly_map, sequence_dir)
    logging.info(f"Proteins info: {proteins_info}")
    database_folder = os.path.join(assembly_folder, "databases")
    if not os.path.isdir(database_folder):
        subprocess.Popen(" ".join(["mkdir ", database_folder]), shell=True)
    os.makedirs(database_folder, exist_ok=True)
    # returns a dictionary with the group number and assemblies in the group
    samples_in_cluster = gc.generate_clusters(
        data,
        args.db_size,
        proteins_info,
        args.study or args.input_dir,
        database_folder,
        sample_assembly_map,
    )
    logging.info(f"Samples in the cluster: {samples_in_cluster}")
    if args.study:
        fd.build_db(
            args.study, database_folder, assembly_folder, samples, samples_in_cluster
        )
    elif args.input_dir:
        fd.build_db(
            "study", database_folder, assembly_folder, samples, samples_in_cluster
        )
    for file in os.listdir(database_folder):
        if file.endswith(".faa") and not file.startswith("unique"):
            protein_file = os.path.join(database_folder, file)
            fd.uniq_proteins(database_folder, file)
            fd.remove_file(protein_file)
    logging.info("Completed")
    logging.info("Runtime is {} seconds".format(time.time() - starttime))


if __name__ == "__main__":
    log_file = "logs/db_generate.log"
    logging.basicConfig(
        filename=log_file,
        filemode="a",
        level=logging.DEBUG,
        format="%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )
    main()
