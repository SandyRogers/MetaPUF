import argparse
import sys
import os
import pandas as pd
import subprocess


def get_args():
    """Collect the inputs for the executable scripts."""
    parser = argparse.ArgumentParser(
        description="""This script is going to map the samples, raw files and searching databases
            between PRIDE and MGnify, then exeucte shell scripts. """
    )

    parser.add_argument('-info', '--sampleinfo', dest='info_file',
                        help="the sample information file, where stores the mapping information for samples, databases and raw files")
    parser.add_argument('-out', '--output', dest='output_folder', help="the folder path for saving files")
    parser.add_argument('-exe', '--exefile', dest='exe_file', help="the exec file")

    args = parser.parse_args()

    if args.info_file is None:
        sys.exit('Error: no info file (sample information file)')
    if args.output_folder is None:
        sys.exit('Error: no input folder path provided!')
    if args.exe_file is None:
        sys.exit('Error: no exec file provided!')

    return args

# sample_info = pd.read_csv("sample_info_final.csv", sep=',')
def thermorawfileparser(exe_file, info_file, output_folder):
    """
    downloading raw file via FTP and convert them into peak list files based on the sample information
    :exe_file:      the thermo raw file parser exeucte file (whole path)
    :info_file:     sample information file contains raw files, samples, mapped searching databases
    :output_folder: input/output folder path for raw/peak list data
    """
    sample_info = pd.read_csv(info_file, sep=',')
    sample_info = sample_info[['Sample','Raw file','Raw file URLs']].drop_duplicates()
    RawURLs  = sample_info['Raw file URLs'].dropna().to_list()

    if len(RawURLs) > 0:
        samples = {k: g["Raw file URLs"].tolist() for k,g in sample_info.groupby("Sample")}
        for sample in samples.keys():
            # generate the folder if it is not there
            folder = os.path.join(output_folder, sample)
            commandline = "mkdir -p " + folder + ";"
            # download the raw files if they are not provided locally.
            for url in samples[sample]:
                commandline += " wget -P " + folder + " -i " + url + ";"

            commandline += " mono " + exe_file + " -d=" + folder + " -o=" + folder
            # -f=1 means mzML format file
            commandline += " -f=0 -m=0 &> logs/" + sample + "_thermorawfileparser.log"
            subprocess.run(commandline, shell=True)

    else:
        samples = {k: g["Raw file"].tolist() for k,g in sample_info.groupby("Sample")}
        for sample in samples.keys():
            # generate the folder if it is not there
            folder = os.path.join(output_folder, sample)
            commandline = "mkdir -p " + folder + ";"
            commandline += " mono " + exe_file + " -d=" + folder + " -o=" + folder
            # -f=1 means mzML format file
            commandline += " -f=0 -m=0 &> logs/" + sample + "_thermorawfileparser.log"
            subprocess.run(commandline, shell=True)


def main():
    """
    Main function
    """
    args = get_args()

    thermorawfileparser(args.exe_file, args.info_file, args.output_folder)


if __name__ == "__main__":
    main()
