import argparse
import sys
import pandas as pd
import subprocess


def get_args():
    """Collect the inputs for the executable scripts."""
    parser = argparse.ArgumentParser(
        description="""This script is going to map the samples, raw files and searching databases
            between PRIDE and MGnify, then exeucte shell scripts. """
    )

    parser.add_argument('-in', '--input', dest='input_file',
                        help="the sample information file, where stores the mapping information for samples, databases and raw files")
    parser.add_argument('-out', '--output', dest='output_folder', help="the folder path for saving output files")
    parser.add_argument('-jar', '--javafile', dest='jar_file', help="the jar file")
    parser.add_argument('-par', '--parameters', dest='par_file', help="the parameter file for searchgui script")
    parser.add_argument('-fn', '--firstname', dest='first_name', help="contact first name")
    parser.add_argument('-ln', '--lastname', dest='last_name', help="contact last name")
    parser.add_argument('-ce', '--contactemail', dest='email', help="contact email")
    parser.add_argument('-ca', '--contactaddress', dest='address', help="contact address")
    parser.add_argument('-on', '--orgname', dest='org_name', help="organization name")
    parser.add_argument('-oe', '--orgemail', dest='org_email', help="organization email")
    parser.add_argument('-oa', '--orgaddress', dest='org_address', help="organization address")
    parser.add_argument('-s', '--searchgui', dest='searchgui',
                        help="seargui search process", action='store_true', default=False)
    parser.add_argument('-p', '--peptideshaker', dest='peptideshaker',
                        help="peptideshaker load process", action='store_true', default=False)
    args = parser.parse_args()

    if args.input_file is None:
        sys.exit('Error: no input file (sample information file)')
    if args.output_folder is None:
        sys.exit('Error: no output path provided!')
    if args.jar_file is None:
        sys.exit('Error: no jar file provided!')
    if not (args.searchgui or args.peptideshaker):
        sys.exit('Error: no processing statues is given: searchgui or peptideshaker')
    if args.searchgui + args.peptideshaker > 1:
        sys.exit('Error: more than one processing statues is provided')

    return args


def searchgui_search(jar_file, input_file, output_folder, par_file):
    """
    run searchgui shell script within a loop based on the mapping sample information
    :input_file:    sample information file contains raw files, samples, mapped searching databases
    :jar_file:      the searchgui jar file
    :output_folder: output folder path for searchgui search results
    :par_file:      the searchgui parameters file
    """
    sample_info = pd.read_csv(input_file, sep=',')
    samples = sample_info['Sample'].drop_duplicates().to_list()

    for sample in samples:
        params = " -spectrum_files"
        params += " input/Raw/" + sample
        params += " -fasta_file"
        fasta = sample_info.loc[sample_info['Sample'] == sample, 'Db_name'].iloc[0]
        fasta = fasta.split('.')[0]
        params += " assemblies/databases/" + fasta + "_concatenated_target_decoy.fasta"
        params += " -output_folder"
        params += " " + output_folder
        params += " -id_params"
        params += " " + par_file
        params += " -xtandem 1 -msgf 1"
        params += " -output_default_name"
        params += " " + sample + "_searchgui"
        params += " -output_option 0 -output_data 1 -output_date 0"
        params += " -log logs/SearchGUI_search &> "
        params += "logs/" + sample + "_SearchGUI_search.log"
        params += " && touch "
        params += "searchgui/" + sample + "_searchgui.zip"

        commandline = "java -cp " + jar_file + " eu.isas.searchgui.cmd.SearchCLI"
        commandline += params

        # print(commandline)
        subprocess.run(commandline, shell=True)


def peptideshaker_load(jar_file, input_file, output_folder, first_name, last_name,
                        email, address, org_name, org_email, org_address):
    """
    run peptideshaker shell script within a loop based on the mapping sample information
    :input_file:    sample information file contains raw files, samples, mapped searching databases
    :jar_file:      the peptideshaker jar file
    :output_folder: output folder path for peptideshaker results
    """
    sample_info = pd.read_csv(input_file, sep=',')
    samples = sample_info['Sample'].drop_duplicates().to_list()

    for sample in samples:
        # params = " -reference " + sample
        params = " -experiment 'peptideshaker' -sample " + sample + " -replicate 1 "
        params += " -identification_files searchgui/" + sample + "_searchgui.zip"
        params += " -out_reports " + output_folder
        params += " -reports 6,9"
        # params += " -report_prefix " + sample + "_"
        params += " -output_file peptideshaker/" + sample + "_peptideshaker.mzid"
        params += " -contact_first_name '" + first_name + "'"
        params += " -contact_last_name '" + last_name + "'"
        params += " -contact_email '" + email + "'"
        params += " -contact_address '" + address + "'"
        params += " -organization_name '" + org_name + "'"
        params += " -organization_email '" + org_email + "'"
        params += " -organization_address '" + org_address + "'"
        params += " &> logs/" + sample + "_PeptideShaker_load.log"

        commandline = "java -cp " + jar_file + " eu.isas.peptideshaker.cmd.PeptideShakerCLI"
        commandline += params

        subprocess.run(commandline, shell=True)


def main():
    """
    Main function
    """
    args = get_args()

    if args.searchgui:
        searchgui_search(args.jar_file, args.input_file, args.output_folder, args.par_file)
    elif args.peptideshaker:
        peptideshaker_load(args.jar_file, args.input_file, args.output_folder, args.first_name, args.last_name,
                            args.email, args.address, args.org_name, args.org_email, args.org_address)
    else:
        sys.exit("BUG! this should not happen.")


if __name__ == "__main__":
    main()
