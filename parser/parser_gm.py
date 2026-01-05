#!/usr/bin/env python3
# author: Laroussi Bachri
# December 29th , 2025
# M1 BBS , University of Toulouse
# Bioinformatics project for genomics
import re
import sys

"""
Module to parse GeneMark output and generate a GFF3 file.
The GFF3 format is described here:
https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
We used this resource to build the GFF3 file.
"""

def extract_cds_genemark(input_file : str):
    """
    Docstring for extract_cds_genemark
    Extracts CDS information from a GeneMark file.
    Args:
        input_file (str): Path to the GeneMark input file.
    Returns:
        list: A list of lists, each sub-list containing the fields of a CDS.
    """

    cds_list = []
    state = "outside" # initial state: outside the CDS section, which section we are in

    with open(input_file, "r") as f:
        for line_traitee in f:
            line = line_traitee.strip()  # remove spaces at the beginning/end of the line

            if state == "outside":  # looking for the CDS header
                if "List of Open reading frames" in line:
                    state = "header" # found
                continue

            if state == "header":
                # waiting for the "bar" under the column header
                if line.startswith("--------"):
                    state = "data"
                continue


            if "List of Regions of interest" in line or "ABOUT THE MATRIX USED" in line: # end of CDS section
                break # exit loop; no useful information afterward

            if line == "" or line.startswith("Left") or line.startswith("end") or line.startswith("--------"): # ignore empty or uninformative lines
                continue 

            # CDS lines start with a number (after strip) and contain direct/complement + fr
            if line[0].isdigit() and (("direct" in line) or ("complement" in line)) and ("fr" in line):
                fields = line.split()  # split on spaces by default, obtain as list
                # fields attendu: [left, right, strand, 'fr', frame, coding_prob, start_prob]
                # Ex: ['2402','3313','direct','fr','2','0.47','0.06']
                cds_list.append(fields)
            
    return cds_list

def extract_info_Genemark(input_file : str):
    """
    Docstring for extract_info_Genemark
    Extracts the sequence identifier, source program name, and sequence length from a GeneMark output file.

    Args:
        input_file (str): Path to the GeneMark input file.
    return:

    seqid :  Sequence ID (e.g., 1_5404 or FASTA fragment name)

    source :  Program name (GeneMark, GeneMark.hmm, scan_for_matches)
    """
    seqid_pat = r'Sequence:\s(\S*)' # pattern to extract sequence ID
    taille_pat = r'Sequence length:\s(\d+)' # pattern to extract sequence length
    id = ''
    source = 'Unknown' 
    taille = '.'
    with open(input_file , 'r') as fh:
        for line in fh:
            if line.startswith('Sequence:'):
                seqid_match = re.search(seqid_pat , line) 
                if seqid_match:
                    id = seqid_match.group(1) # extract ID
                else:
                    id = 'unknown_id'
            if 'GENEMARK' in line.upper(): # detect GeneMark program
                # The file does not explicitly show the name of the program used
                source = 'GeneMark' # assume GeneMark if we find this word
            if line.startswith('Sequence length:'):
                taille_match = re.search(taille_pat , line)
                if taille_match:
                    taille = int(taille_match.group(1))
            if id != '' and source != 'Unknown' and taille != '.':
                break # found both required pieces of information
    return id, source , taille
            
def write_gff3(input_file : str, output_file : str):
    """
    Docstring for write_gff3
    This function writes a GFF3 file from data extracted from a GeneMark file.
    It uses the `extract_info_Genemark` and `extract_cds_genemark` functions to obtain the required information.
    
    :param input_file: path to the input GeneMark file
    :param output_file: path to the output GFF3 file
    :return: confirmation message for writing the .GFF3
    """

    with open(output_file, 'w') as out_fh:
        # Write GFF3 header
        seq_id , source , taille = extract_info_Genemark(input_file)
        cds_list = extract_cds_genemark(input_file)
        out_fh.write("##gff-version 3\n")
        out_fh.write(f'##sequence-region {seq_id} 1 {taille}\n')
        # Write each CDS in GFF3 format
        for idx, fields in enumerate(cds_list):
            start = min(int(fields[0]), int(fields[1]))
            end = max(int(fields[0]), int(fields[1]))
            feature_type = "CDS"
            strand = '+' if fields[2] == "direct" else '-'
            phase = str(int(fields[4]) - 1)  # fr 1->0, fr 2->1, fr3->2
            score = fields[5] if len(fields) > 5 else '.'
            attributes = f"ID=CDS_{idx+1};Note={source}_prediction"

            gff_line = f"{seq_id}\t{source}\t{feature_type}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{attributes}\n"
            out_fh.write(gff_line)
    return f"GFF3 file '{output_file}' written successfully."

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: ./parser_gm.py <input_file.gm> <output_file.gff3>")
        sys.exit(1)
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    write_gff3(input_file , output_file)
