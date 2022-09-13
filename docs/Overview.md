# **Overview**

This workflow offers the users to integrate publicly available metagenomics, metatranscriptomics and metaproteomics datasets on PRIDE and MGnify portals.



Software |  Version |  Purpose
--|---|--
  Snakemake| 7.3.8  |  Wrap the pipelines
Sourmash  | 4.4.2  |  Create DNA sketches to group similar assemblies together
 mg-toolkit |0.10.1   |  Pull information from MGify web portal
 ThermoRawFileParser | 1.2.3 | Convert Thermo Raw files into different formats of spectral files.
 Mono | 6.12.0 | Wrapper around the .net (C#) ThermoFisher ThermoRawFileReader library for running on Linux
 SearchGUI | 4.0.41 | Combine multiple search engines to perform peptide identification from a protein sequence database.
 PeptideShaker | 2.0.33 | Interpret proteomics identification results from multiple search and generate peptide and protein reports.



The workflow is defined in snakemake and is available at [link to github](https://github.com/PRIDE-reanalysis/MetaPUF.git)
