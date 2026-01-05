#!/usr/bin/env python3
# author: Laroussi Bachri
# December 29th , 2025
# M1 BBS , University of Toulouse
# Bioinformatics project for genomics
import re
import sys

"""
Module to parse GeneMark.hmm output and generate a GFF3 file.
The GFF3 format is described here:
https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
We used this resource to build the GFF3 file.
"""
def extract_infos_GenemarkHMM(input_file: str):
    """
    Docstring for extract_infos_GenemarkHMM
    Extracts the sequence identifier, program name, and CDS information
    from a genemarkhmm output file.
    :param input_file: path to the genemarkhmm output file
    :return: tuple (seq_id, source, cds_list)
        seq_id : sequence ID (e.g., FASTA fragment name)
        source : program name (e.g., GeneMark.hmm)
        cds_list : list of CDS entries extracted from the file
    """
    seq_id = ''
    source = ''
    cds_list = []
    with open(input_file , 'r') as fh:
        state = 'outside' # initial state before the data section
        for line_traitee in fh:
            line = line_traitee.strip()
            if 'Version' in line: # line containing the program name
                source = line.split()[0]
            if line.startswith('FASTA definition line:'): # line containing the sequence ID
                seq_match = re.search(r':\s(\S*)' , line)
                if seq_match:
                    seq_id = seq_match.group(1)
            if line.startswith('Predicted genes'): # start of the data section
                state = 'data'
                continue

            if state == 'data' and line and line[0].isdigit(): # verify line starts with a digit
                fields = line.split() # extract fields
                cds_list.append(fields) 
    
    return seq_id , source , cds_list

def taille_seq(seq_fasta:str):
    """
    Docstring for taille_seq
    
    :param seq_fasta: fasta file
    :return: sequence length
    """
    taille = 0
    with open(seq_fasta , 'r') as f:
        for line in f:
            if not line.startswith(">"):
                taille += len(line.strip())
    return taille



def nettoyer_int(x: str):
    """
    Docstring for nettoyer_int
    Avoid conversion errors when turning a string into an integer by removing >/<.
    (Notably on positions in genemarkhmm)
    :param x: string
    return: int
    """
    return int(x.replace("<","").replace(">",""))



def write_gff3(input_file: str , output_file: str ,seq_fasta=None):
    """
    Docstring for write_gff3
    Write a GFF3 file from genemarkhmm output
    :param input_file: genemarkhmm output file
    :param output_file: GFF3 output file
    :param seq_fasta: optional fasta file to obtain the sequence length if available
    
    This function uses extract_infos_GenemarkHMM and taille_seq to obtain the necessary information.

    return: confirmation message
    """
    seq_id , source , cds_list= extract_infos_GenemarkHMM(input_file)
    taille = None
    if seq_fasta: # optional
        taille = taille_seq(seq_fasta)

    with open(output_file, 'w') as out_fh:
        out_fh.write("##gff-version 3\n")

        if taille: # optional if length known
            out_fh.write(f'##sequence-region {seq_id} 1 {taille}\n') # indicate sequence length
        
        for i , fields in enumerate(cds_list):
            left = nettoyer_int(fields[2]) 
            right = nettoyer_int(fields[3])
            start = min(left, right) # start <= end mandatory in GFF3 according to documentation
            end = max(left, right)

            score = '.' # no score in genemarkhmm
            strand = fields[1] # + or -

            if strand == '+': # forward strand
                phase = (start - 1) % 3 # calculate phase from start position
            else:
                phase = (end - 1) % 3 # reverse strand; calculate from end position

            classe = fields[5] if len(fields)>=6 else '.' # between typical and atypical (horizontal transfer, for example)
            attributs = f"ID=CDS_{i+1};Note={source}_prediction;Class={classe}"

            out_fh.write(f"{seq_id}\t{source}\tCDS\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{attributs}\n")
    return f"GFF3 file '{output_file}' written successfully."



if __name__ == '__main__':
    
    if len(sys.argv) < 3:
        print("Usage: ./parser_genemarkhmm.py <input_file_gmhmmp.out> <output_file.gff3> [seq_fasta]")
        sys.exit(1)
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    fasta = sys.argv[3] if len(sys.argv) ==4 else None 
    if len(sys.argv) == 4: # optionnel
        seq_fasta = sys.argv[3]
        write_gff3(input_file , output_file , seq_fasta)
    else:   
        write_gff3(input_file , output_file)

    
