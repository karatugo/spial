# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 19:21:41 2020

@author: bselcuk
"""

def conservation_check(place_num,align_dict):
    place_num=int(place_num)
    place_count=0
    residue_count=0
    result_aa_list=[]
    result_aa_number=[]
    total=0
    count=0
    for seq in align_dict:
        if count==0:
            reference_seq=align_dict[seq]
            while place_count<place_num:
                amino_acid=reference_seq[place_count]
                if amino_acid!="-":
                    residue_count+=1
                place_count+=1
            count=1 
        total+=1
        target_aa=align_dict[seq][place_count-1]
        if target_aa=="-":
            total-=1
        if target_aa not in result_aa_list:
            result_aa_list.append(target_aa)
            result_aa_number.append(0)
        idx=result_aa_list.index(target_aa)
        result_aa_number[idx]+=1
    max_num=max(result_aa_number)
    idx_max=result_aa_number.index(max_num)
    aa=result_aa_list[idx_max]
    if total==0:
        return [0,"-",-1,"none"]
    #print ("The conservation for {} th residue: {:.2f}%".format(residue_num+1,100*accuracy/total))
    return [max_num*100/total,aa,residue_count]

def alignment_dict(MA_file):
    result_dict={}
    alignment_file=open(MA_file,"r")
    while 1:
        line=alignment_file.readline()
        if line=="":
            break
        if line[0]==">":
            header=line.strip()
            result_dict[header]=""
        else:
            result_dict[header]+=line.strip()
    length=len(result_dict[header])
    return (result_dict,length)

def spial(alignment,groupA,groupB,specificity_threshold,structral_or_sequence):
    result_dict={}
    alignment_d=alignment_dict(alignment)
    align_length=alignment_d[1]
    alignmentA_dict=fasta_divider(groupA,alignment_d[0])
    alignmentB_dict=fasta_divider(groupB,alignment_d[0])
    for i in range(align_length):
        result_A=conservation_check(i+1,alignmentA_dict)
        conservationA=result_A[0]
        aa_A=result_A[1]
        residue_num_A=result_A[2]
        result_B=conservation_check(i+1,alignmentB_dict)
        conservationB=result_B[0]
        aa_B=result_B[1]
        if structral_or_sequence=="structural":
            place=GPCR_convert(keys[i],gene)
        elif structral_or_sequence=="sequence":
            place=residue_num_A
        if aa_A=="-":
            continue
        if conservationA>=specificity_threshold and conservationB<specificity_threshold:
            if aa_A!=aa_B:
                if conservationA not in result_dict:
                    result_dict[conservationA]=[]
                result_dict[conservationA].append([place,aa_A,i+1])
        elif conservationA>=specificity_threshold and aa_A!=aa_B:
            if conservationA not in result_dict:
                result_dict[conservationA]=[]
            result_dict[conservationA].append([place,aa_A,i+1])
    return result_dict

def spial_result_parse(result_dict,structural_or_sequence):
    conservation_list=list(result_dict.keys())
    conservation_list.sort(reverse=True)
    PyMOLinput=""
    for value in conservation_list:
        for result in result_dict[value]: 
            place_in_alignment=result[0]
            PyMOLinput+=str(place_in_alignment)+"+"
            print ("Conservation:{:.2f}\tResidue number:{}\tAmino acid:{}\tMSA Position:{}".format(value,result[0],result[1],result[2]))
    print ("-----------------------------")
    return PyMOLinput[:-1]

def GPCR_convert(res,gene):
    import pickle
    gene=gene.lower()
    with open(r"D:\Users\suuser\Desktop\{}_notation_dict".format(gene), 'rb') as handle:
        conv_dict = pickle.load(handle)
        if res not in conv_dict:
            return "-"
    return conv_dict[res]


def fasta_divider(gene_list,align_dict):
    result_dict={}
    for gene in gene_list:
        count=0
        for seq in align_dict:
            if seq.split("|")[-2]=="{}_HUMAN".format(gene):
                result_dict[seq]=align_dict[seq]
                count=1
                continue
            if count==1:
                if seq.split("|")[-1]=="9606":
                    break
                else:
                    result_dict[seq]=align_dict[seq]
    return result_dict

gene="ADRB1" #The protein that you are finding its specific residues
keys=['1x25', '1x26', '1x27', '1x28', '1x29', '1x30', '1x31', '1x32', '1x33', '1x34', '1x35', '1x36', '1x37', '1x38', '1x39', '1x40', '1x41', '1x42', '1x43', '1x44', '1x45', '1x46', '1x47', '1x48', '1x49', '1x50', '1x51', '1x52', '1x53', '1x54', '1x55', '1x56', '1x57', '1x58', '1x59', '1x60', '12x48', '12x49', '12x50', '12x51', '2x37', '2x38', '2x39', '2x40', '2x41', '2x42', '2x43', '2x44', '2x45', '2x46', '2x47', '2x48', '2x49', '2x50', '2x51', '2x52', '2x53', '2x54', '2x55', '2x551', '2x56', '2x57', '2x58', '2x59', '2x60', '2x61', '2x62', '2x63', '2x64', '2x65', '2x66', '23x49', '23x50', '23x51', '23x52', '3x21', '3x22', '3x23', '3x24', '3x25', '3x26', '3x27', '3x28', '3x29', '3x30', '3x31', '3x32', '3x33', '3x34', '3x35', '3x36', '3x37', '3x38', '3x39', '3x40', '3x41', '3x42', '3x43', '3x44', '3x45', '3x46', '3x47', '3x48', '3x49', '3x50', '3x51', '3x52', '3x53', '3x54', '3x55', '3x56', '34x50', '34x51', '34x52', '34x53', '34x54', '34x55', '34x56', '34x57', '4x38', '4x39', '4x40', '4x41', '4x42', '4x43', '4x44', '4x45', '4x46', '4x47', '4x48', '4x49', '4x50', '4x51', '4x52', '4x53', '4x54', '4x55', '4x56', '4x57', '4x58', '4x59', '4x60', '4x61', '4x62', '4x63', '4x64', '45x50', '45x51', '45x52', '5x36', '5x37', '5x38', '5x39', '5x40', '5x41', '5x42', '5x43', '5x44', '5x45', '5x46', '5x461', '5x47', '5x48', '5x49', '5x50', '5x51', '5x52', '5x53', '5x54', '5x55', '5x56', '5x57', '5x58', '5x59', '5x60', '5x61', '5x62', '5x63', '5x64', '5x65', '5x66', '5x67', '5x68', '5x69', '5x70', '5x71', '5x72', '5x73', '5x74', '6x24', '6x25', '6x26', '6x27', '6x28', '6x29', '6x30', '6x31', '6x32', '6x33', '6x34', '6x35', '6x36', '6x37', '6x38', '6x39', '6x40', '6x41', '6x42', '6x43', '6x44', '6x45', '6x46', '6x47', '6x48', '6x49', '6x50', '6x51', '6x52', '6x53', '6x54', '6x55', '6x56', '6x57', '6x58', '6x59', '6x60', '6x61', '7x30', '7x31', '7x32', '7x33', '7x34', '7x35', '7x36', '7x37', '7x38', '7x39', '7x40', '7x41', '7x42', '7x43', '7x45', '7x46', '7x47', '7x48', '7x49', '7x50', '7x51', '7x52', '7x53', '7x54', '7x55', '8x47', '8x48', '8x49', '8x50', '8x51', '8x52', '8x53', '8x54', '8x55', '8x56', '8x57', '8x58', '8x59', '5x75', '5x76', '2x67', '7x56', '5x411', '6x22', '6x23', '8x60', '8x61', '8x62', '8x63', '8x64', '8x65', '8x66', '8x67', '8x68', '8x69', '6x21', '4x65', '4x66', '4x67', '4x68', '4x69', '5x32', '5x33', '5x34', '5x35', '7x27', '7x28', '7x29', '1x22', '1x23', '1x24', '5x31']
path=r"D:\Users\suuser\Desktop\Structral_alignments\\"
#print (spial(path+"ADRB2_structuralMA.fasta",path+"ADRB1-3_structuralMA.fasta",95))
structure_result=spial(path+"ADRB1-3_structuralMA.fasta",["ADRB1"],["ADRB2","ADRB2"],95,"structural")
sequence_result=spial(r"D:\Users\suuser\Desktop\ADRB1-3_alligned.fasta",["ADRB1"],["ADRB2","ADRB2"],95,"sequence")
print ("PyMOL input:",spial_result_parse(structure_result,"structural"))
print ("====================================")
print ("PyMOL input:",spial_result_parse(sequence_result,"sequence")+"\n")