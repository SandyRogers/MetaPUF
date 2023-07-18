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


import os
import pandas as pd

pd.options.mode.chained_assignment = None  # default='warn'


def calculate_max_count(df: pd.DataFrame) -> int:
    """Generate maximum spectrum count value from a dataframe
    :param pd.DataFrame df: input dataframe
    :return int: maximum spectruem count value
    """
    sc_values = []
    for i in range(len(df)):
        if not ";" in str(df["Spectrum Counting"][i]):
            sc_values.append(df["Spectrum Counting"][i])
    max_spectrum_count = max(sc_values)
    return max_spectrum_count


def protein_report_processing(
    protein_report, protein_list: list, pride_id: str
) -> pd.DataFrame:
    """
    processing of peptide reports to yield a dataframe
    :param pride_id: ID of the associated metaproteomics study
    :param protein_list: list of expressed proteins
    :param protein_report: processed protein report

    :return pd.DataFrame: returns two dataframes
    """

    # create two empty dataframe
    pep_info_high = pd.DataFrame(
        columns=[
            "digest",
            "contig_name",
            "protein_start",
            "protein_end",
            "strand",
            "Attributes",
        ]
    )
    # pep_info_low=pd.DataFrame(columns=['digest','contig_name','protein_start','protein_end','strand','Attributes'])
    for item in protein_list:
        item = str(item)
        all_info = []
        temp_data = (protein_report.loc[protein_report["digest"] == item]).reset_index(
            drop=True
        )
        unambiguous_peptides = []
        ambiguous_peptides = []
        sc = []
        # low_confidence_ambiguous_peptides=[]
        # low_confidence_all_info=[]
        # Unique_peptide_to_protein_mapping="None"
        # Ambiguous_peptide_to_protein_mapping="None"
        for i in range(len(temp_data)):
            # I have comment this because we don't need to check whether the protein belongs to one protein groups or not,
            # we only check if the peptide is unique mapped to the protein, which means we only check the #Proteins.
            if int(temp_data["#Proteins"][i]) == 1:
                unambiguous_peptides.append(
                    temp_data["Sequence"][i]
                    + " [PSMs:"
                    + str(temp_data["Validated PSMs"][i])
                    + "]"
                )

            else:
                ambiguous_peptides.append(
                    temp_data["Sequence"][i]
                    + " [PSMs:"
                    + str(temp_data["Validated PSMs"][i])
                    + "]"
                )

            sc.append(temp_data["Spectrum Counting"][i])
            contig = "_".join((temp_data["contig_name"][i]).split("_")[:-1])
            start = temp_data["protein_start"][i]
            end = temp_data["protein_end"][i]
            strand = temp_data["strand"][i]

        if len(set(unambiguous_peptides)) == 0 and len(set(ambiguous_peptides)) >= 1:
            all_info.append(
                (
                    str(item),
                    contig,
                    start,
                    end,
                    strand,
                    "ID="
                    + str(item)
                    + ";type=Protein;Unique_peptide_to_protein_mapping=None;Ambiguous_peptide_to_protein_mapping=True;ambiguous_sequences="
                    + ",".join(ambiguous_peptides)
                    + ";pride_id="
                    + pride_id
                    + ";semiquantitative_expression_spectrum_count="
                    + str(sc),
                )
            )
        elif len(set(unambiguous_peptides)) >= 1 and len(set(ambiguous_peptides)) == 0:
            all_info.append(
                (
                    str(item),
                    contig,
                    start,
                    end,
                    strand,
                    "ID="
                    + str(item)
                    + ";type=Protein;Unique_peptide_to_protein_mapping=True;unambiguous_sequences="
                    + ",".join(unambiguous_peptides)
                    + ";Ambiguous_peptide_to_protein_mapping=None;pride_id="
                    + pride_id
                    + ";semiquantitative_expression_spectrum_count="
                    + str(sc),
                )
            )
        else:
            all_info.append(
                (
                    str(item),
                    contig,
                    start,
                    end,
                    strand,
                    "ID="
                    + str(item)
                    + ";type=Protein;Unique_peptide_to_protein_mapping=True;unambiguous_sequences="
                    + ",".join(unambiguous_peptides)
                    + ";Ambiguous_peptide_to_protein_mapping=True;ambiguous_sequences="
                    + ",".join(ambiguous_peptides)
                    + ";pride_id="
                    + pride_id
                    + ";semiquantitative_expression_spectrum_count="
                    + str(sc),
                )
            )

        df_high = pd.DataFrame(
            all_info,
            columns=[
                "digest",
                "contig_name",
                "protein_start",
                "protein_end",
                "strand",
                "Attributes",
            ],
        )
        pep_info_high = pd.concat([pep_info_high, df_high], ignore_index=True)

        # get the max spectrum value for the processed protein file
        max_spectrum_count = protein_report["max_spectrum_count"].unique()

    return pep_info_high, max_spectrum_count


def gff_generation_unique(
    attributes_file: str, assembly_name: str, out_folder: str, max_spectrum_count
):
    """
    function to generate the gff file from the data file
    :param str attributes_file: input dataframe with all attributes for the gff file
    :param str assembly_name: name of the assemble from the study
    :param str out_folder: filepath of the output folder
    :param _type_ max_spectrum_count: value of max spectral count in the given assembly

    """
    gff_data = attributes_file[
        ["contig_name", "protein_start", "protein_end", "strand", "Attributes"]
    ]
    gff_data["strand"] = gff_data["strand"].map({-1: "-", 1: "+"})
    gff_data = gff_data.rename(
        columns={
            "contig_name": "seqid",
            "protein_start": "start",
            "protein_end": "end",
            "strand": "strand",
            "Attributes": "attributes",
        }
    )
    gff_data = gff_data.assign(source="PeptideShaker")
    gff_data = gff_data.assign(type="CDS")
    gff_data = gff_data.assign(score=".")
    gff_data = gff_data.assign(phase=".")
    cols = gff_data.columns.to_list()
    col = (
        cols[:1] + cols[5:7] + cols[1:3] + cols[7:8] + cols[3:4] + cols[-1:] + cols[4:5]
    )
    gff_data = gff_data[col]
    topline = "##gff-version 3"
    df_flatten = list(zip(*map(gff_data.get, gff_data)))
    out_file = os.path.join(out_folder, assembly_name + "_unique_peptides.gff")
    with open(out_file, "w", buffering=1) as out_handle:
        print("##gff-version 3", file=out_handle)
        print(
            "##max_spectrum_count_value_in_study=" + str(max_spectrum_count[0]),
            file=out_handle,
        )
        for row in df_flatten:
            print("\t".join([str(val) for val in row]), file=out_handle)


def gff_generation_ambiguous(attributes_file: str, assembly_name: str, out_folder: str):
    gff_data = attributes_file[
        ["contig_name", "protein_start", "protein_end", "strand", "Attributes"]
    ]
    gff_data["strand"] = gff_data["strand"].map({-1: "-", 1: "+"})
    gff_data = gff_data.rename(
        columns={
            "contig_name": "seqid",
            "protein_start": "start",
            "protein_end": "end",
            "strand": "strand",
            "Attributes": "attributes",
        }
    )
    gff_data = gff_data.assign(source="PeptideShaker")
    gff_data = gff_data.assign(type="CDS")
    gff_data = gff_data.assign(score=".")
    gff_data = gff_data.assign(phase=".")
    cols = gff_data.columns.to_list()
    col = (
        cols[:1] + cols[5:7] + cols[1:3] + cols[7:8] + cols[3:4] + cols[-1:] + cols[4:5]
    )
    gff_data = gff_data[col]
    topline = "##gff-version 3"
    df_flatten = list(zip(*map(gff_data.get, gff_data)))
    out_file = os.path.join(out_folder, assembly_name + "_ambiguous_peptides.gff")
    with open(out_file, "w", buffering=1) as out_handle:
        print("##gff-version 3", file=out_handle)
        for row in df_flatten:
            print("\t".join([str(val) for val in row]), file=out_handle)
