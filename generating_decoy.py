import argparse
import sys
import os
import pandas as pd
import subprocess


def get_args():
    """Collect the inputs for the executable scripts."""
    parser = argparse.ArgumentParser(
        description="""This script is going to add the human and contanminat protein sequences
                into the generated databases and also add decoy into every database"""
    )

    parser.add_argument('-info', '--sampleinfo', dest='info_file',
                        help="the final sample information file, where also stores generaed databases")
    parser.add_argument('-human', '--humandb', dest='human_db', help="the human reference protein sequences")
    parser.add_argument('-crap', '--crapdb', dest='crap_db', help="the contanminat sequences")
    parser.add_argument('-jar', '--jarfile', dest='jar_file', help="the jar file")
    parser.add_argument('-tmp', '--tmpdir', dest='tmp_dir', help="the temprare folder storing the flag file for the process")

    args = parser.parse_args()

    if args.info_file is None:
        sys.exit('Error: no info file (sample information file)')
    if args.human_db is None:
        sys.exit('Error: no human reference protein database provided!')
    if args.crap_db is None:
        sys.exit('Error: no contanminat protein database provided!')
    if args.exe_file is None:
        sys.exit('Error: no jar file provided!')
    if args.tmp_dir is None:
        sys.exit('Error: no tmp folder provided!')

    return args


def searchgui_decoy(jar_file, info_file, human_db, crap_db, tmpdir):
    """
    adding human, contanminat and decoy to the protein databases.
    :jar_file:      the searchgui jar file (whole path)
    :info_file:     final sample information file also contains mapped searching databases
    :human_db:      the human protein reference sequences
    :crap_db:       the contanminat protein sequences
    """
    sample_info = pd.read_csv(info_file, sep=',')
    proteins = sample_info['Db_name'].drop_duplicates().to_list()
    filenames = [human_db, crap_db]

    for protein in proteins:
        protein_file = " assemblies/databases/" + protein
        with open(protein_file, 'w') as outfile:
            for fname in filenames:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)

        params = " -in" + protein_file
        params += " -decoy &> "
        params += "logs/" + protein + "_SearchGUI_search.log"

        commandline += "java -cp " + jar_file + " eu.isas.searchgui.cmd.FastaCLI"
        commandline += params

        subprocess.run(commandline, shell=True)

    commandline = "mkdir -p " + tmpdir
    subprocess.run(commandline, shell=True)
    flag_txt = "proteins_decoy_generated_check.txt"
    with open(flag_txt, 'w') as f:
        f.write('All protein databases have been processed!')


def main():
    """
    Main function
    """
    args = get_args()

    searchgui_decoy(args.jar_file, args.info_file, args.human_db, args.crap_db, args.tmpdir)


if __name__ == "__main__":
    main()
