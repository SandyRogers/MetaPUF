import argparse
import json
import logging
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd


def get_args():
    """Collect the inputs for the executable scripts."""
    parser = argparse.ArgumentParser(
        description="""This script is going to add the human and contanminat protein sequences
                into the generated databases and also add decoy into every database"""
    )

    parser.add_argument(
        "-info",
        "--sampleinfo",
        dest="info_file",
        help="the final sample information file, where also stores generaed databases",
    )
    parser.add_argument(
        "-human",
        "--humandb",
        dest="human_db",
        help="the human reference protein sequences",
    )
    parser.add_argument(
        "-crap", "--crapdb", dest="crap_db", help="the contanminat sequences"
    )
    parser.add_argument("-jar", "--jarfile", dest="jar_file", help="the jar file")
    parser.add_argument(
        "-sgpars",
        "--searchgui_params",
        dest="searchgui_params",
        help="json-like string of options to pass to searchgui",
    )
    parser.add_argument(
        "-adir",
        "--assemblies_dir",
        dest="assemblies_dir",
        help="the output assemblies directory path",
    )

    args = parser.parse_args()

    if args.info_file is None:
        sys.exit("Error: no info file (sample information file)")
    if args.human_db is None:
        sys.exit("Error: no human reference protein database provided!")
    if args.crap_db is None:
        sys.exit("Error: no contanminat protein database provided!")
    if args.jar_file is None:
        sys.exit("Error: no jar file provided!")

    return args


def searchgui_decoy(
    jar_file, info_file, human_db, crap_db, searchgui_params, assemblies_folder
):
    """
    adding human, contaminant and decoy to the protein databases.
    :jar_file:      the searchgui jar file (whole path)
    :info_file:     final sample information file also contains mapped searching databases
    :human_db:      the human protein reference sequences
    :crap_db:       the contanminat protein sequences
    :assemblies_folder: the output/assemblies dir path
    """
    sample_info = pd.read_csv(info_file, sep=",")
    proteins = sample_info["Db_name"].dropna().drop_duplicates().to_list()

    config = json.loads(searchgui_params)
    SEARCHGUI_PAR_PARAMS = " ".join(
        [
            "-%s %s" % (k, "'%s'" % v if isinstance(v, str) else str(v))
            for k, v in config.items()
        ]
    )

    for protein in proteins:
        protein_file = f"{assemblies_folder}/databases/{protein}"
        print(f"On {protein = } with {protein_file = }")
        fasta = Path(protein_file).resolve().with_suffix(".fasta").as_posix()
        par = f"{Path(protein_file).parent.resolve()}/{Path(protein_file).stem}_searchgui.par"
        db = f"{Path(fasta).parent.resolve()}/{Path(fasta).stem}_concatenated_target_decoy.fasta"
        log = Path(protein_file).stem + "_SearchGUI_params.log"

        commandline = (
            "cat " + protein_file + " " + human_db + " " + crap_db + " > " + fasta + ";"
        )
        # generating the protein decoy
        params = " -in " + fasta
        params += " -decoy &> "
        params += "logs/" + protein.split(".")[0] + "_SearchGUI_search.log;"

        commandline += " java -cp " + jar_file + " eu.isas.searchgui.cmd.FastaCLI"
        commandline += params

        # generating the protein parameter file
        params = f" -out {par} "
        params += f" -db {db} "
        params += SEARCHGUI_PAR_PARAMS
        params += f" &> logs/{log};"
        commandline += (
            " java -cp "
            + jar_file
            + " eu.isas.searchgui.cmd.IdentificationParametersCLI"
        )
        commandline += params

        print(f"Will run {commandline}")
        subprocess.run(commandline, shell=True)
        # print(commandline)

    # commandline = "mkdir -p " + s_dir
    # subprocess.run(commandline, shell=True)
    flag_txt = (
        f"{assemblies_folder}/databases/proteins_decoy_params_generated_check.txt"
    )
    with open(flag_txt, "w") as f:
        f.write("All protein databases have been processed!")


def main():
    """
    Main function
    """
    args = get_args()

    searchgui_decoy(
        args.jar_file,
        args.info_file,
        args.human_db,
        args.crap_db,
        args.searchgui_params,
        args.assemblies_dir,
    )


if __name__ == "__main__":
    main()
