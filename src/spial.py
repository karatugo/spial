import sys
import os
import json
import re
import operator
from .fasta_dict import *

base_path = "./uploaded-fa-files/"
fasta_files_path = base_path 

def get_aa_conservation(fasta_file):
    aa_count_dict = {}
    fasta_dict = get_fasta_dict(base_path,fasta_file)
    length = len(fasta_dict[list(fasta_dict.keys())[0]])
    #print(length)
    i = 0
    for pos in range(0,length):
        aa_count_dict[i] = {}
        for k,v in fasta_dict.items():
            if not aa_count_dict:
                aa_count_dict[i].update({v[pos]:0})
            if v[pos] not in aa_count_dict[i].keys():
                #print(v[pos])
                aa_count_dict[i].update({v[pos]:0})
                #print(aa_count_dict)
            if v[pos] in  aa_count_dict[i].keys():
                aa_count_dict[i][v[pos]] += 1
              
        #print(i)
        i += 1
    #print(aa_count_dict)
    conservation_dict = {}
    for position, counts in aa_count_dict.items():
        conservation_dict[position] = {}
        raw_length = float(len(fasta_dict.keys()))
        max_score_aa = max(aa_count_dict[position].items(), key=operator.itemgetter(1))[0]
        max_score = float(aa_count_dict[position][max_score_aa]/raw_length)
        conservation_dict[position].update({max_score_aa:max_score})
    #print(conservation_dict)
    return conservation_dict

def get_sdp(conservation_dict1,conservation_dict2,consensus_treshold,specificity_treshold):
    sdp_dict = {}
    consensus_treshold = float(consensus_treshold)
    specificity_treshold = float(specificity_treshold)
    #print(len(list(conservation_dict1.keys())))
    for position in range(0,len(list(conservation_dict1.keys()))):
        sdp_dict[position] = {}
        aa1 = list(conservation_dict1[position].keys())[0]
        aa2 = list(conservation_dict2[position].keys())[0]
        score1 = conservation_dict1[position][aa1]
        score2 = conservation_dict2[position][aa2]
       
        if aa1 == aa2:
            if score1 > consensus_treshold and score2 > consensus_treshold:
                sdp_dict[position] = {aa1:{'score':score1,'type':'consensus'}}
            elif score1 < consensus_treshold and score2 < consensus_treshold:
                sdp_dict[position] = {aa1:{'score':score1,'type':'non-specific'}}
            elif score1 > consensus_treshold and score2 < consensus_treshold:
                if score1 > specificity_treshold:
                    sdp_dict[position] = {aa1:{'score':score1,'type':'specific'}}
                else:
                    sdp_dict[position] = {aa1:{'score':score1,'type':'non-specific'}}
        elif aa1 != aa2:
            if score1 > specificity_treshold and score2 > specificity_treshold:
                sdp_dict[position] = {aa1:{'score':score1,'type':'separately specific'}}
            elif score1 > specificity_treshold and score2 < specificity_treshold:
                sdp_dict[position] = {aa1:{'score':score1,'type':'specific'}}

           
    #print(sdp_dict)
    return sdp_dict

def dump_to_json(sdp_dict,json_file):
    with open(base_path +json_file, 'w') as f:
        json.dump(sdp_dict, f, ensure_ascii=False, indent=4)
   


if __name__ == "__main__":
    fasta_file1= sys.argv[1] 
    fasta_file2 = sys.argv[2] 
    consensus_treshold = sys.argv[3]
    specificity_treshold = sys.argv[4]
    conservation_dict1 = get_aa_conservation(fasta_file1)
    conservation_dict2 = get_aa_conservation(fasta_file2)
    sdp_dict = get_sdp(conservation_dict1,conservation_dict2,consensus_treshold,specificity_treshold)
    dump_to_json(sdp_dict,"sdp.json")