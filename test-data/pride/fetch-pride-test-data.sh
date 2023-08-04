#!/bin/bash
cd $(dirname "$0")
mkdir -p S6
curl -o S6/S6.raw ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2017/07/PXD005780/S6.raw
mkdir -p S13
curl -o S13/S13.raw ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2017/07/PXD005780/S13.raw
