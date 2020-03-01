import os
import sys
def get_fasta_dict(path,fasta_file):
    fasta_dict = {}
    
    with open(path + fasta_file,'r') as f:
        for line in f:
            if line.startswith('>'):
                header = line.strip()
                fasta_dict[header] = ''
            else:
                fasta_dict[header] += line.strip()
            
    f.close()

    return fasta_dict
    