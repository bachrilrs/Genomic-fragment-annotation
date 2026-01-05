#!/usr/bin/env python3
# author: Laroussi Bachri
# December 29th , 2025
# M1 BBS , University of Toulouse
# Bioinformatics project for genomics

import re
import sys
"""
Module to parse Scan_for_matches output and generate a GFF3 file.
The GFF3 format is described here:
https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
We used this resource to build the GFF3 file.
"""
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

def parse_scanformatches(input_file: str):
    """
    Docstring for parse_scanformatches
    Extracts match information from a scan_for_matches file.
    Args:
        input_file (str): Path to the scan_for_matches input file.
    Returns:
        list: A list of lists, each sub-list containing the fields of a match.
    """
    seq_id = None
    positions = {}
    i=0
    with open(input_file, "r") as f:
        for line_traitee in f:
            line = line_traitee.strip()  # remove spaces at the beginning/end of the line

            if line.startswith(">") and seq_id is None:
                seq_id_match = re.search(r'^>(\S+):' , line)
                seq_id = seq_id_match.group(1)

            if line.startswith(">"):
                i+=1
                pos_match = re.search(r':\[(\d*),(\d*)]', line)
                if not pos_match:
                    continue
                positions[i] = [pos_match.group(1),pos_match.group(2)]
            start_match = re.search(r'\b(atg|ttg|gtg)\b' , line)
            if start_match:
                positions[i].append(start_match.group(1))

    return seq_id , positions

def write_gff3(input_file: str , output_file: str , feature_type , seq_fasta=None):
    """
    Docstring for write_gff3
    Write a GFF3 file from scan_for_matches output
    :param input_file: scan_for_matches output file
    :param output_file: GFF3 output file
    :param seq_fasta: optional fasta file to obtain the sequence length if available

    return: confirmation message
    """
    seq_id , positions = parse_scanformatches(input_file)
    source = 'scan_for_matches'
    taille = None
    feature_type = feature_type  # or RBS/promoter/terminator depending on context
    if seq_fasta: # optional
        taille = taille_seq(seq_fasta)

    with open(output_file, 'w') as out_fh:
        out_fh.write("##gff-version 3\n")

        if taille: # optional if length known
            out_fh.write(f'##sequence-region {seq_id} 1 {taille}\n') # indicate sequence length
        
        for i  in sorted(positions):
            codon_start = None
            if len(positions[i]) == 3:
                start , end, codon_start= positions[i]
            else:
                start , end= positions[i]


            gff_start = min(int(start), int(end))
            gff_end = max(int(start), int(end)) # start <= end mandatory in GFF3 according to documentation

            score = '.' # no score in scan_for_matches
            strand = '+' if int(start) < int(end) else '-'

            phase = '.' # no phase for RBS/promoter/terminator
            if codon_start:
                attributs = f"ID={feature_type}_{i};Note={source}_prediction;Start_Codon={codon_start}"
            else:
                attributs = f"ID={feature_type}_{i};Note={source}_prediction"
            out_fh.write(f"{seq_id}\t{source}\t{feature_type}\t{gff_start}\t{gff_end}\t{score}\t{strand}\t{phase}\t{attributs}\n")
    return f"GFF3 file '{output_file}' written successfully."

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Usage: ./parser_scanformatches.py <input_file> <output_file.gff3> <feature_type> [seq_fasta]")
        sys.exit(1)
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    feature_type = sys.argv[3]  # RBS/promoter/terminator
    fasta = sys.argv[4] if len(sys.argv) == 5 else None 
    if len(sys.argv) == 5: # optionnel
        write_gff3(input_file , output_file, feature_type , fasta)
    else:   
        write_gff3(input_file , output_file,feature_type)
