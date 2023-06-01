# MetaPUF
____________________________________

####  An approach to integrate metagenomics, metatranscriptomics and metaproteomics data found in public resources such as MGnify (for metagenomics/metatranscriptomics) and the PRIDE database (for metaproteomics). When these omics techniques are applied to the same sample, their integration offers new opportunities to understand the structure (metagenome) and functional expression (metatranscriptome and metaproteome) of the microbiome.

# Installation
____________________________________

Please follow the instructions from the [tutorial](https://metapuf-tutorial.readthedocs.io/en/latest/installation.html#installation) 

The packages and their versions:
- Python3
- Snakemake==7.3.8
- biopython==1.79
- mg-toolkit==0.10.1
- pandas==1.4.2
- sourmash==4.4.2

The versions of the tools ( [`ThermoRawFileParser`](https://github.com/compomics/ThermoRawFileParser), `SearchGui` and `PeptideShaker` ) that we used in our pipeline are not updated to the latest ones, you can still use the latest version of ThermoRawFileParser, but you will need to download other two tools from these links: [PeptideShaker 1.16.45](https://genesis.ugent.be/maven2/eu/isas/peptideshaker/PeptideShaker/) and [SearchGui 3.3.20](https://genesis.ugent.be/maven2/eu/isas/searchgui/SearchGUI/), after downloading, please unzip them and move the folders under the path of `workflow/bin`.
Please find more detailed documentation from [MetaPUF](https://metapuf-tutorial.readthedocs.io/en/latest/index.html)

## (Linux) Requirements
[Mono](https://www.mono-project.com/download/stable/#download-lin) (install mono-complete if you encounter "assembly not found" errors).


# Example Usage
____________________________________

- The example folder has a sample of post processed reports metaproteomics study (marine?) and sample_info.txt that contains all the information for generating metagenomics assemblies and build the protein sequence databases.

- The user can update the size of the protein serach database in the config.proteomics.yaml file in the config folder or else it takes the default size of 1GB.

- The user can update the config.proteomics.yaml file in the  config.proteomics.yaml file in the config folder with the output results folder and the name of the PRIDE ID.

- Since the pipeline is a long-running task, it is recommended to use some terminal multiplexer such as `screen` or other job control tools to run the pipeline in the background, and the memory for running the pipeline should be big enough as well. 

## Tips for running Snakemake
- You can run a dry-run to check for any syntax errors 
```
 $ Snakemake  -np
```

- To run the workflow
```
 $ Snakemake --cores 4
```


- Tips: IF the pipeline got collapsed during running, you can always try to run a dry-run `Snakemake  -np` first to check how many rules have been successful executed, and if you are sure that some files are generated correctly, you can use `snakemake --cleanup-metadata <filenames>` to skip these files to be re-generated. However, sometimes `snakemake --cleanup-metadata <filenames>` doesn't work, you can also try to manually delete the `.snakemake/incomplete` directory.


# Core contributors and collaborators

# Code of Conduct
As part of our efforts toward delivering open and inclusive science, we follow the [Contributor Covenant Code of Conduct for Open Source Projects](https://www.contributor-covenant.org/version/2/0/code_of_conduct/).

# How to cite

# Copyright notice
