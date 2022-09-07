# -*- coding: utf-8 -*-

# Copyright 2022 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import os
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'


def calculate_max_count(reports_dir):
    sc_values=[]
    for file_in in os.listdir(reports_dir):
        if file_in.endswith('_peptide_report.csv'):
            file_name=os.path.join(reports_dir, file_in)
            sample = (pd.read_csv(file_name, sep=',',index_col=[0])).reset_index()
            for i in range(len(sample)):
                if not ";" in str(sample['Spectrum Counting'][i]):
                    sc_values.append(sample['Spectrum Counting'][i])
    max_spectrum_count=max(sc_values)
    return max_spectrum_count


def protein_report_processing(protein_report, protein_list: list,pride_id: str):
    """
    processing of peptide reports to yield a dataframe
    :param pride_id: ID of the associated metaproteomics study
    :param protein_list: list of expressed proteins
    :param protein_report: processed protein report
    """
    all_info=[]
    pep_info_high=pd.DataFrame(columns=['digest','contig_name','protein_start','protein_end','strand','Attributes'])
    pep_info_low=pd.DataFrame(columns=['digest','contig_name','protein_start','protein_end','strand','Attributes'])
    for item in protein_list:
        item = str(item)
        all_info=[]
        temp_data=(protein_report.loc[protein_report['digest']==item]).reset_index(drop=True)
        unambiguous_peptides=[]
        ambiguous_peptides=[]
        low_confidence_ambiguous_peptides=[]
        low_confidence_all_info=[]
        Unique_peptide_to_protein_mapping="None"
        Ambiguous_peptide_to_protein_mapping="None"
        for i in range(len(temp_data)):
            if (int(temp_data["#Proteins"][i] == 1) and int(temp_data["Validated Protein Groups"][i]==1)):
                unambiguous_peptides.append(temp_data["Sequence"][i]+" [PSMs:"+str(temp_data["Validated PSMs"][i])+"]")
                sc=temp_data["Spectrum Counting"][i]
                contig="_".join((temp_data["contig_name"][i]).split("_")[:-1])
                start=temp_data["protein_start"][i]
                end=temp_data["protein_end"][i]
                strand=temp_data["strand"][i]
            elif (int(temp_data["#Proteins"][i] > 1) and int(temp_data["Validated Protein Groups"][i]==1)):
                ambiguous_peptides.append(temp_data["Sequence"][i]+" [PSMs:"+str(temp_data["Validated PSMs"][i])+"]")
                sc=temp_data["Spectrum Counting"][i]
                contig="_".join((temp_data["contig_name"][i]).split("_")[:-1])
                start=temp_data["protein_start"][i]
                end=temp_data["protein_end"][i]
                strand=temp_data["strand"][i]
            else:
                low_confidence_ambiguous_peptides.append(temp_data["Sequence"][i]+" [PSMs:"+str(temp_data["Validated PSMs"][i])+"]")
                contig="_".join((temp_data["contig_name"][i]).split("_")[:-1])
                start=temp_data["protein_start"][i]
                end=temp_data["protein_end"][i]
                strand=temp_data["strand"][i]
        unambiguous_peptides = list(set(unambiguous_peptides))
        ambiguous_peptides = list(set(ambiguous_peptides))
        low_confidence_ambiguous_peptides=list(set(low_confidence_ambiguous_peptides))
        if len(unambiguous_peptides)==0 and len(ambiguous_peptides)>=1:
            all_info.append( (str(item) ,contig,start,end,strand, "ID="+str(item)+";type=Protein;Unique_peptide_to_protein_mapping=None;Ambiguous_peptide_to_protein_mapping=True;ambiguous_sequences="+",".join(ambiguous_peptides)+";pride_id="+pride_id+";semiquantitative_expression_spectrum_count="+str(sc)) )
        elif len(unambiguous_peptides)>=1 and len(ambiguous_peptides)==0:
            all_info.append( (str(item) ,contig,start,end,strand, "ID="+str(item)+";type=Protein;Unique_peptide_to_protein_mapping=True;unambiguous_sequences="+",".join(unambiguous_peptides)+";Ambiguous_peptide_to_protein_mapping=None;pride_id="+pride_id+";semiquantitative_expression_spectrum_count="+str(sc)) )
        else:
            all_info.append( (str(item) ,contig,start,end,strand, "ID="+str(item)+";type=Protein;Unique_peptide_to_protein_mapping=True;unambiguous_sequences="+",".join(unambiguous_peptides)+";Ambiguous_peptide_to_protein_mapping=True;ambiguous_sequences="+",".join(ambiguous_peptides)+";pride_id="+pride_id+";semiquantitative_expression_spectrum_count="+str(sc)) )
        df_high = pd.DataFrame(all_info, columns=['digest', 'contig_name','protein_start','protein_end','strand','Attributes'])
        pep_info_high = pd.concat([pep_info_high,df_high], ignore_index=True)
        if len(low_confidence_ambiguous_peptides)>=1:
            low_confidence_all_info.append( (str(item) ,contig,start,end,strand, "ID="+str(item)+";type=Protein;Unique_peptide_to_protein_mapping=None;Ambiguous_peptide_to_protein_mapping=True;ambiguous_sequences="+",".join(ambiguous_peptides)+";pride_id="+pride_id) )
        df_low = pd.DataFrame(low_confidence_all_info, columns=['digest', 'contig_name','protein_start','protein_end','strand','Attributes' ])
        pep_info_low = pd.concat([pep_info_low,df_low], ignore_index=True)
    return pep_info_high,pep_info_low

def gff_generation_high(report_dir:str, attributes_file: str, assembly_name:str, out_folder: str):
    max_spectral_value = calculate_max_count(report_dir)
    gff_data=attributes_file[['contig_name','protein_start','protein_end','strand','Attributes']]
    gff_data['strand'] = gff_data['strand'].map({-1: '-', 1: '+'})
    gff_data = gff_data.rename(columns={'contig_name':'seqid','protein_start':'start','protein_end':'end','strand':'strand','Attributes':'attributes'})
    gff_data = gff_data.assign(source='PeptideShaker')
    gff_data = gff_data.assign(type='CDS')
    gff_data = gff_data.assign(score='.')
    gff_data = gff_data.assign(phase='.')
    cols = gff_data.columns.to_list()
    col=cols[:1] + cols[5:7] +cols[1:3] +cols[7:8] +cols[3:4] +cols[-1:] +cols[4:5]
    gff_data=gff_data[col]
    topline='##gff-version 3'
    df_flatten = list(zip(*map(gff_data.get, gff_data)))
    out_file=os.path.join(out_folder, assembly_name+"_high_confidence.gff")
    with open(out_file,'w', buffering=1) as out_handle:
        print('##gff-version 3', file=out_handle)
        print('##max_sprectrum_count_value_in_study='+str(max_spectral_value), file=out_handle)
        for row in df_flatten:
            print('\t'.join([str(val) for val in row]), file=out_handle)


def gff_generation_low(attributes_file: str, assembly_name:str, out_folder: str):
    gff_data=attributes_file[['contig_name','protein_start','protein_end','strand','Attributes']]
    gff_data['strand'] = gff_data['strand'].map({-1: '-', 1: '+'})
    gff_data = gff_data.rename(columns={'contig_name':'seqid','protein_start':'start','protein_end':'end','strand':'strand','Attributes':'attributes'})
    gff_data = gff_data.assign(source='PeptideShaker')
    gff_data = gff_data.assign(type='CDS')
    gff_data = gff_data.assign(score='.')
    gff_data = gff_data.assign(phase='.')
    cols = gff_data.columns.to_list()
    col=cols[:1] + cols[5:7] +cols[1:3] +cols[7:8] +cols[3:4] +cols[-1:] +cols[4:5]
    gff_data=gff_data[col]
    topline='##gff-version 3'
    df_flatten = list(zip(*map(gff_data.get, gff_data)))
    out_file=os.path.join(out_folder, assembly_name+"_low_confidence.gff")
    with open(out_file,'w', buffering=1) as out_handle:
        print('##gff-version 3', file=out_handle)
        for row in df_flatten:
            print('\t'.join([str(val) for val in row]), file=out_handle)
