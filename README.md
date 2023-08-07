# MetaPUF
____________________________________

####  An approach to integrate metagenomics, metatranscriptomics and metaproteomics data found in public resources such as MGnify (for metagenomics/metatranscriptomics) and the PRIDE database (for metaproteomics). When these omics techniques are applied to the same sample, their integration offers new opportunities to understand the structure (metagenome) and functional expression (metatranscriptome and metaproteome) of the microbiome.

# Installation
You need a working [installation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) of [Snakemake](https://snakemake.readthedocs.io/en/stable/).
Then:
```shell
git clone <this-repo>
```

The pipeline uses conda environments to manage dependencies, which are handled automatically if you run snakemake with the `--use-conda` flag.

It also relies on some tools ([`ThermoRawFileParser`](https://github.com/compomics/ThermoRawFileParser), [`SearchGui`](http://compomics.github.io/projects/searchgui) and [`PeptideShaker`](http://compomics.github.io/projects/peptide-shaker))
which do not have conda packages or docker images available for the versions we used.
These tools are downloaded on-the-fly by snakemake, so you do not need to install them separately.

# Example usage (test data-set)
There is a small test-data set, using a few assemblies from MGnify and two RAW files from PRIDE.
To fetch the (~GB size) RAW files, which are too big for this git repository:
```shell
./test-data/pride/fetch-pride-test-data.sh
```
This downloads two RAW files into `test-data/pride/`.

Then:
```shell
conda activate snakemake # (assuming you installed snakemake with conda, into an env called snakemake)
cd MetaPUF
snakemake --cores 4 --use-conda
```
This will run the pipeline on the small dataset, and put results into `../test-run`.

# Real usage (configuration)
Edit the `config/config.proteomics.yaml` and `sample_info.csv` files to point the pipeline at real data.
`sample_info.csv` is the mapping of MGnify to PRIDE datasets, and in the config `parameters.input_dir` and `parameters.raw_dir` refer to the MGnify and PRIDE data folders respectively.


## Tips for running Snakemake
- You can run a dry-run to check for any syntax errors
```
 Snakemake  -np
```

- To run the workflow
```
 Snakemake --cores 4 --use-conda
```

- Using LSF on an HPC cluster:
```shell
bsub -n 4 -R "rusage[mem=4096]" -J metapuf -u $USER -o job.log -e job.err snakemake --cores 4 --use-conda
```

- Tips: IF the pipeline got collapsed during running, you can always try to run a dry-run `Snakemake  -np` first to check how many rules have been successful executed, and if you are sure that some files are generated correctly, you can use `snakemake --cleanup-metadata <filenames>` to skip these files to be re-generated. However, sometimes `snakemake --cleanup-metadata <filenames>` doesn't work, you can also try to manually delete the `.snakemake/incomplete` directory.

# Development installation
```shell
git clone https://github.com/PRIDE-reanalysis/MetaPUF.git
cd MetaPUF
conda activate snakemake  # or another conda/venv if you prefer
pip install ".[dev,docs]"
pre-commit install
```

This installs the development requirements, and installs the pre-commit hooks which format the code correctly while commiting changes.
You can also manually format the code using `black .`.
It also installs `mkdocs`, which is used to build the documentation.

## Editing docs
Change the markdown files in the `docs/` folder.
Then `mkdocs serve` to view the documentation site locally.

# Core contributors and collaborators

# Code of Conduct
As part of our efforts toward delivering open and inclusive science, we follow the [Contributor Covenant Code of Conduct for Open Source Projects](https://www.contributor-covenant.org/version/2/0/code_of_conduct/).

# How to cite

# Copyright notice
