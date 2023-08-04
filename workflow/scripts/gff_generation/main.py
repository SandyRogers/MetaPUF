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
import time
from argparse import ArgumentParser

import gff_builder as gb
import pandas as pd


def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)


def main():  # noqa: C901
    """
    Aggregate information from metaproteomics and metagenomics about the expressed proteins
    and generate GFF format files for visualisation
    """
    parser = ArgumentParser(
        description="Aggregate information from metaproteomics and metagenomics about the expressed proteins and generate GFF format files for visualisation"
    )
    parser.add_argument(
        "-i",
        "--sample_info",
        type=str,
        required=True,
        help="Absolute path of the sample metadata file",
    )
    parser.add_argument(
        "-r",
        "--reports_dir",
        type=dir_path,
        required=True,
        help="full path of the directory containing processed metaproteomics reports",
    )
    parser.add_argument(
        "-m",
        "--metag_dir",
        type=dir_path,
        required=True,
        help="full path of the directory containing contig metadata file",
    )
    parser.add_argument(
        "-p",
        "--pride_id",
        type=str,
        required=True,
        help="PRIDE id for the reanalysed study",
    )

    starttime = time.time()
    args = parser.parse_args()
    # read in sample_info file and get the sample names in a list
    sample_info = pd.read_csv(args.sample_info, sep=",")

    # check if the results folder exists or create the results folder
    results_folder = os.path.join(args.reports_dir, "results")
    os.makedirs(results_folder, exist_ok=True)

    # get the filenames as list and concatenate the expressed peptides from all the samples
    sample_file_list = [
        args.reports_dir + "/" + f
        for f in os.listdir(args.reports_dir)
        if f.endswith("_peptide_report.csv")
    ]
    csv_list = []

    for file in sorted(sample_file_list):
        csv_list.append(pd.read_csv(file))

    csv_merged = pd.concat(csv_list, ignore_index=True)
    csv_merged.to_csv(os.path.join(results_folder, "allsamples.csv"), index=False)
    proteins = list(set(csv_merged["Protein"].to_list()))
    # create a new dataframe with unique peptide digests
    unique_proteins = pd.DataFrame(proteins, columns=["digest"])
    # get the list of distinct assemblies in the study
    assemblies = list(set(sample_info["Assembly"].to_list()))

    for assembly in assemblies:
        temp_assembly = pd.read_csv(
            os.path.join(args.metag_dir, assembly + "_contig_info.txt"), sep="\t"
        )
        temp_assembly.columns = [
            "assembly",
            "contig_name",
            "digest",
            "partial_info",
            "protein_start",
            "protein_end",
            "strand",
            "caller",
        ]
        # get all the matching expressed proteins from the contigs_info.txt
        assembly_subset_of_proteins = pd.merge(
            unique_proteins, temp_assembly, on="digest", how="inner"
        )

        # get all the peptides expressed in the given assembly
        assembly_expressed_proteins = assembly_subset_of_proteins.merge(
            csv_merged, left_on="digest", right_on="Protein", how="inner"
        )

        if len(assembly_expressed_proteins) > 0:
            assembly_expressed_proteins["max_spectrum_count"] = gb.calculate_max_count(
                assembly_expressed_proteins
            )
            assembly_expressed_proteins.to_csv(
                os.path.join(results_folder, assembly + "_expressed_proteins.csv")
            )
            expressed_proteins = list(set(assembly_expressed_proteins["digest"]))
            attributes_file, max_spectrum_count = gb.protein_report_processing(
                assembly_expressed_proteins, expressed_proteins, args.pride_id
            )
            if len(attributes_file) >= 1:
                gb.gff_generation_unique(
                    attributes_file, assembly, results_folder, max_spectrum_count
                )

    logging.info("Completed")
    logging.info("Runtime is {} seconds".format(time.time() - starttime))


if __name__ == "__main__":
    log_file = "logs/gff_generate.log"
    logging.basicConfig(
        filename=log_file,
        filemode="a",
        level=logging.DEBUG,
        format="%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )
    main()
