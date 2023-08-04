import argparse
import subprocess
import sys

import pandas as pd


def get_args():
    """Collect the inputs for the executable scripts."""
    parser = argparse.ArgumentParser(
        description="""This script is going to map the samples, raw files and searching databases
            between PRIDE and MGnify, then execute shell scripts. """
    )

    parser.add_argument(
        "-in",
        "--input",
        dest="input_file",
        help="the sample information file, where stores the mapping information for samples, databases and raw files",
    )
    parser.add_argument(
        "-raw",
        "--raw_folder",
        dest="raw_folder",
        help="the folder path where .raw files are found",
        required=False,
    )
    parser.add_argument(
        "-out",
        "--output_folder",
        dest="output_folder",
        help="the folder path in which searchgui and peptideshaker folder will be saved/read",
    )
    parser.add_argument("-jar", "--javafile", dest="jar_file", help="the jar file")
    parser.add_argument(
        "-s",
        "--searchgui",
        dest="searchgui",
        help="seargui search process",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-p",
        "--peptideshaker",
        dest="peptideshaker",
        help="peptideshaker load process",
        action="store_true",
        default=False,
    )
    args = parser.parse_args()

    if args.input_file is None:
        sys.exit("Error: no input file (sample information file)")
    if args.output_folder is None:
        sys.exit("Error: no output path provided!")
    if args.jar_file is None:
        sys.exit("Error: no jar file provided!")
    if not (args.searchgui or args.peptideshaker):
        sys.exit("Error: no processing statues is given: searchgui or peptideshaker")
    if args.searchgui + args.peptideshaker > 1:
        sys.exit("Error: more than one processing statues is provided")

    return args


def searchgui_search(jar_file, input_file, raw_folder, output_folder):
    """
    run searchgui shell script within a loop based on the mapping sample information
    :input_file:    sample information file contains raw files, samples, mapped searching databases
    :jar_file:      the searchgui jar file
    :raw_folder:    folder containing the .raw PRIDE files
    :output_folder: output folder path for searchgui search results
    """
    sample_info = pd.read_csv(input_file, sep=",")
    samples = sample_info["Sample"].drop_duplicates().to_list()

    for sample in samples:
        params = " -spectrum_files"
        params += f" {raw_folder}/" + sample
        # params += " -fasta_file"
        fasta = sample_info.loc[sample_info["Sample"] == sample, "Db_name"].iloc[0]
        print(f"On sample {sample = } with {fasta =}")
        fasta = fasta.split(".")[0]
        params += f" -output_folder {output_folder}/searchgui"
        params += " -id_params"
        params += f" {output_folder}/assemblies/databases/" + fasta + "_searchgui.par"
        params += " -xtandem 1 -msgf 1"
        params += " -output_default_name"
        params += " " + sample + "_searchgui"
        params += " -output_option 0 -output_data 1 -output_date 0"
        params += " -log logs/SearchGUI_search &> "
        params += "logs/" + sample + "_SearchGUI_search.log"
        params += " && touch "
        params += f"{output_folder}/searchgui/{sample}_searchgui.zip"

        commandline = "java -cp " + jar_file + " eu.isas.searchgui.cmd.SearchCLI"
        commandline += params

        print(f"Will run: {commandline}")
        subprocess.run(commandline, shell=True)


def peptideshaker_load(
    jar_file,
    input_file,
    output_folder,
):
    """
    run peptideshaker shell script within a loop based on the mapping sample information
    :input_file:    sample information file contains raw files, samples, mapped searching databases
    :jar_file:      the peptideshaker jar file
    :output_folder: output folder path for peptideshaker results
    """
    sample_info = pd.read_csv(input_file, sep=",")
    samples = sample_info["Sample"].drop_duplicates().to_list()

    for sample in samples:
        print(f"Load peptideshaker for sample {sample}")
        params = " -experiment 'peptideshaker' -sample " + sample + " -replicate 1 "
        params += (
            f" -identification_files {output_folder}/searchgui/{sample}_searchgui.zip"
        )
        params += f" -out_reports {output_folder}/peptideshaker"
        params += " -reports 6,9"
        params += (
            f" -output_file {output_folder}/peptideshaker/{sample}_peptideshaker.mzid"
        )
        params += " &> logs/" + sample + "_PeptideShaker_load.log"

        commandline = (
            "java -cp " + jar_file + " eu.isas.peptideshaker.cmd.PeptideShakerCLI"
        )
        commandline += params
        print(f"Will run peptideshaker command {commandline}")

        subprocess.run(commandline, shell=True)


def main():
    """
    Main function
    """
    args = get_args()

    if args.searchgui:
        searchgui_search(
            args.jar_file, args.input_file, args.raw_folder, args.output_folder
        )
    elif args.peptideshaker:
        peptideshaker_load(
            args.jar_file,
            args.input_file,
            args.output_folder,
        )
    else:
        sys.exit("BUG! this should not happen.")


if __name__ == "__main__":
    main()
