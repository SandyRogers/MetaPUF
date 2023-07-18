import argparse
import sys
import os
import pandas as pd
import subprocess
import yaml


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
        "-p", "--params", dest="params_file", help="the parameters for SearchGUI"
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
    if args.params_file is None:
        sys.exit("Error: no parameters provided!")

    return args


# info_file = 'assemblies/sample_info_final.csv'
# human_db = 'resources/human_db.fa'
# crap_db = 'resources/crap_db.fa'
def searchgui_decoy(jar_file, info_file, human_db, crap_db, params_file):
    """
    adding human, contanminat and decoy to the protein databases.
    :jar_file:      the searchgui jar file (whole path)
    :info_file:     final sample information file also contains mapped searching databases
    :human_db:      the human protein reference sequences
    :crap_db:       the contanminat protein sequences
    """
    sample_info = pd.read_csv(info_file, sep=",")
    proteins = sample_info["Db_name"].drop_duplicates().to_list()

    SEARCHGUI_PAR_PARAMS = ""
    with open(params_file, "r") as yamlfile:
        config = yaml.load(yamlfile, Loader=yaml.FullLoader)
        SEARCHGUI_PAR_PARAMS = " ".join(
            [
                "-%s %s" % (k, "'%s'" % v if isinstance(v, str) else str(v))
                for k, v in config["searchgui"]["par"].items()
            ]
        )
        yamlfile.close()

    for protein in proteins:
        protein_file = "assemblies/databases/" + protein
        fasta = protein_file.split(".")[0] + ".fasta"
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
        params = " -out " + protein_file.split(".")[0] + "_searchgui.par"
        params += " -db " + fasta.split(".")[0] + "_concatenated_target_decoy.fasta "
        params += SEARCHGUI_PAR_PARAMS
        params += " &> logs/" + protein.split(".")[0] + "_SearchGUI_params.log;"
        commandline += (
            " java -cp "
            + jar_file
            + " eu.isas.searchgui.cmd.IdentificationParametersCLI"
        )
        commandline += params

        subprocess.run(commandline, shell=True)
        # print(commandline)

    # commandline = "mkdir -p " + s_dir
    # subprocess.run(commandline, shell=True)
    flag_txt = "assemblies/databases/proteins_decoy_params_generated_check.txt"
    with open(flag_txt, "w") as f:
        f.write("All protein databases have been processed!")


def main():
    """
    Main function
    """
    args = get_args()

    searchgui_decoy(
        args.jar_file, args.info_file, args.human_db, args.crap_db, args.params_file
    )


if __name__ == "__main__":
    main()
