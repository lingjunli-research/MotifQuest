# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 14:41:13 2024

@author: lafields2
"""

import pandas as pd
import csv
from Bio.SeqIO.FastaIO import SimpleFastaParser

partial_motif_path = r"D:\Manuscripts\2024_MotifAutomation\Algorithm_Validation\Databases\AMP\partial_motif_0.40_AMP_scores_2.csv"
full_motif_path = r"D:\Manuscripts\2024_MotifAutomation\Algorithm_Validation\Databases\AMP\full_motif_0amp_scores_2.csv"
fasta_path = r"D:\Manuscripts\2024_MotifAutomation\Algorithm_Validation\Databases\AMP\AMP_database_UNMC.fasta"
output_path = r"D:\Manuscripts\2024_MotifAutomation\Algorithm_Validation\Databases\AMP\processed_dbs"
#manual_db_path = r"D:\Manuscripts\2024_MotifAutomation\Algorithm_Validation\Databases\Manual_NP_motif_DB.csv"
min_len = 2

aa_list = ['A','G','I','L','P','V','F','W','Y','D','E','R','H','K','S','T','C','M','N','Q']

scores_rep = pd.read_csv(partial_motif_path)

fasta_to_df = []
title_to_df = []

with open(fasta_path) as fasta_file:  # Will close handle cleanly
    for title, sequence in SimpleFastaParser(fasta_file):
        fasta_to_df.append(sequence)

fasta_df = pd.DataFrame()
fasta_df['Sequence'] = fasta_to_df
#fasta_df['Title'] = title_to_df

partial_motif_list = scores_rep['motif'].values.tolist()

all_motif_options_storage = []

for a in partial_motif_list:
    for b in aa_list:
        replaced_a = a.replace('X',b)
        res = any(replaced_a in sub for sub in fasta_to_df)
        
        if res == True:
            if replaced_a not in all_motif_options_storage:
                all_motif_options_storage.append(replaced_a)
            else:
                pass
        else:
            pass

full_motif = pd.read_csv(full_motif_path)
full_motif_list = full_motif['motif'].values.tolist()

all_motifs = full_motif_list + all_motif_options_storage

all_motifs_filtered = []
for x in all_motifs:
    if len(x) >= min_len:
        all_motifs_filtered.append(x)
    else:
        pass

motif_rep_out = pd.DataFrame()
motif_rep_out['Sequence'] = all_motifs_filtered

out = output_path + '\\combined_AMP_motif_DB.csv'
with open(out,'w',newline='') as filec:
        writerc = csv.writer(filec)
        motif_rep_out.to_csv(filec,index=False)
        
partial_motifs_filtered = []
for x in all_motif_options_storage:
    if len(x) >= min_len:
        partial_motifs_filtered.append(x)
    else:
        pass

partial_motif_rep_out = pd.DataFrame()
partial_motif_rep_out['Sequence'] = partial_motifs_filtered

out = output_path + '\\partial_AMP_motif_DB.csv'
with open(out,'w',newline='') as filec:
        writerc = csv.writer(filec)
        partial_motif_rep_out.to_csv(filec,index=False)
        
full_motifs_filtered = []
for x in full_motif_list:
    if len(x) >= min_len:
        full_motifs_filtered.append(x)
    else:
        pass

full_motif_rep_out = pd.DataFrame()
full_motif_rep_out['Sequence'] = full_motifs_filtered

out = output_path + '\\full_AMP_motif_DB.csv'
with open(out,'w',newline='') as filec:
        writerc = csv.writer(filec)
        full_motif_rep_out.to_csv(filec,index=False)

# paired_seq_storage = []
# paired_mot_storage = []



# for value in all_motifs:
#     if len(value) >= min_len:
#         for string in fasta_to_df:
#             if value in string:
#                 paired_seq_storage.append(string)
#                 paired_mot_storage.append(value)

# paired_rep = pd.DataFrame()
# paired_rep['Motif'] = paired_mot_storage
# paired_rep['Sequence'] = paired_seq_storage

# paired_rep_fam = paired_rep.merge(fasta_df, on='Sequence')

# out = output_path + '\\full_NP_motif_DB_3_w_fam_new.csv'
# with open(out,'w',newline='') as filec:
#         writerc = csv.writer(filec)
#         paired_rep_fam.to_csv(filec,index=False)
# # #%%
# manual_db = pd.read_csv(manual_db_path)
# manual_db_motif_list = manual_db['Sequence'].values.tolist()

# man_motif_storage = []
# man_seq_storage = []
# true_storage = []

# for mot in manual_db_motif_list:
#     for seq in fasta_to_df:
#         if mot in seq:
#             man_motif_storage.append(mot)
#             man_seq_storage.append(seq)
        
#         if mot not in seq:
#             pass

# for mot in manual_db_motif_list:       
#     if mot in man_motif_storage:
#         pass
#     elif mot not in man_motif_storage:
#         man_motif_storage.append(mot)
#         man_seq_storage.append('Matching seq not found')     

# for x in man_seq_storage:
#     if x in fasta_to_df:
#         true_storage.append(True)
#     if x not in fasta_to_df:
#         true_storage.append(False)

# man_motif_check = pd.DataFrame()
# man_motif_check['Motif'] = man_motif_storage
# man_motif_check['Sequence'] = man_seq_storage
# man_motif_check['Status'] = true_storage


# out = output_path + '\\manual_motif_db_corrected_new.csv'
# with open(out,'w',newline='') as filec:
#         writerc = csv.writer(filec)
#         man_motif_check.to_csv(filec,index=False)