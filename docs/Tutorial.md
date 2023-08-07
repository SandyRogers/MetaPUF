# Tutorials

Here below you can find tutorials on how to use some of the MetaPUF pipelines and standalones. This gives you a flavor of what can be achieved using MetaPUF.

## Installation

#### Install conda executable

We recommend to install Miniconda Python3 distribution from [Download](https://conda.io/en/latest/miniconda.html).

#### Install Snakemake

We recommend to [install Snakemake using conda](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html), if you don't already have it.

This tutorial assumes you install Snakemake into a conda env called `snakemake`.

#### From GitHub source code

You can install MetaPUF locally from GitHub

```shell
conda activate snakemake
git clone https://github.com/PRIDE-reanalysis/MetaPUF.git
cd MetaPUF
```

The dependencies are handled by Snakemake's creation of conda envs, and by some Snakemake rules that download non-conda dependencies.

## Running a test project

A small test dataset is available consisting of a couple of MGnify analyses and a couple of PRIDE raw files.
The PRIDE raw files are too big for GitHub, so they can be fetched using a script:
```shell
./test-data/pride/fetch-pride-test-data.sh
```
This downloads two RAW files into `test-data/pride/`, which will use around 2GB disk space.

The config file in `config/config.proteomics.yaml` and the sample mapping in `config/sample_info.csv` are already written to work with this dataset.

To run the dataset:

```shell
snakemake --cores 4 --use-conda
```
Outputs will be in `../test-run`, logs will be in `logs/`.

## Inputs

The workflow requires the following inputs in the config.yaml file:

* ` Study: (ERP124921) | input_dir: Absolute path to the input directory`

This parameter takes in European Nucleotide Archive (ENA) secondary study accession: starts with (ERP|DRP|SRP) followed by six digits.

You can also input your own metagenomics and/or metatranscriptomics assemblies and metaproteomics data from the same samples. There are some instructions regarding the naming convention to follow:

- Assembly files should have file names ending with "".
- The workflow accepts protein predicted from Prodigal.
- The protein file names should end with "".

* `Pride_id: PRIDE dataset accession For example (PXD005780)`

Depending on the network, it is better to download the dataset (raw files) beforehand and put them in the correct folder, each sample should have a sub-folder for the Raw files. You can also use the pipeline to download the files for you, but it may take time; also, if the connection failes, all the files will need to be re-downloaded again.

** Example of folder stuctrue: (PXD005780 for example)
Under the root path of the MetaPUF project, create a folder named `input/Raw`, under this path, you need to create sub-folders for each sample: `S06`, `S07`, ... etc. And put the raw files into the correlated sub-folder. If there are multiple raw files for one sample, please put all those raw files into the same sub-folder which represents for one sample.

* `Version: The version of MGnidy analysis.`

This accepts either "4.1" or "5.0" as string values. The default is "5.0".
* `Db_size: Size of the proteins sequence database in bytes.`

The default is 1073741824 (integer data type).
* `outputdir: Directory path to save the output`

The user needs to enter a .csv file that contains all the sample information as shown in the example below

Sample  | PRIDE Accession  | Raw file  | Raw file URLs  | Sample Accession  | Assembly
--|---|---|---|---|--
 S6 |  PXD005780 | S6.raw | https://ftp.pride.ebi.ac.uk/pride/data/archive/2017/07/PXD005780/S6.raw  | ERS1509315  |  ERZ1669330

This workflow looks for this file in the config directory. You don't need to provide the `PRIDE Accession` number and `Raw file URLs` if you have the data locally, and just leave the column as blank, however the header is still needed.

## Quick start

Snakemake dry run command

```
Snakemake -n --cores 4
```

To run the workflow with Snakemake

```
Snakemake --cores 4
```

## Outputs

The outputs include the following files
* `Protein sequence database`

* `Report informing assemblies in each protein sequence db`

* `Protein reports from PeptideShaker`

* `Peptide report from PeptideShaker`

* `Post-processing reports`

* `High confidence GFF format report`
* `Low confidence GFF format report`
