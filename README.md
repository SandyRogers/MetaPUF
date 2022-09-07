# MetaPUF
____________________________________

####  Integration of metagenomics, metatranscriptomics and metaproteomics data in the public domain between the PRIDE database and MGnify
# Installation
____________________________________

The packages and their versions:
- Python3
- Snakemake==7.3.8
- biopython==1.79
- mg-toolkit==0.10.1
- pandas==1.4.2
- sourmash==4.4.2

- install with conda (recommended)

- $ conda create -y -n environment -c conda-forge -f environment.yml

# Example Usage
____________________________________

- The example folder has a sample of post processed reports metaproteomics study (marine?) and sample_info.txt that contains all the information for generating metagenomics assemblies and build the protein sequence databases.

- The user can update the size of the protein serach database in the config.proteomics.yaml file in the config folder or else it takes the default size of 1GB.

- The user can update the config.proteomics.yaml file in the  config.proteomics.yaml file in the config folder with the output results folder and the name of the PRIDE ID.

- You can run a dry-run to check for any syntax errors
- $ Snakemake  -n --cores 4

- To run the workflow
- $ Snakemake --cores 4


# Core contributors and collaborators

# Code of Conduct
As part of our efforts toward delivering open and inclusive science, we follow the [Contributor Covenant Code of Conduct for Open Source Projects](https://www.contributor-covenant.org/version/2/0/code_of_conduct/).

# How to cite

# Copyright notice
