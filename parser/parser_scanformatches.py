#!/usr/bin/env python3
# author: Laroussi Bachri
# December 29th , 2025
# M1 BBS , Université de Toulouse
# Projet Bioinformatique pour la génomique

import re
import sys
"""
Module pour parser la sortie de Scan_for_matches et générer un fichier GFF3.
Le format GFF3 est décrit ici :
https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
Nous avons utilisé cette ressource pour construire le fichier GFF3.
"""
def taille_seq(seq_fasta:str):
    """
    Docstring pour taille_seq
    
    :param seq_fasta: fichier fasta
    :return: taille de la séquence
    """
    taille = 0
    with open(seq_fasta , 'r') as f:
        for line in f:
            if not line.startswith(">"):
                taille += len(line.strip())
    return taille

def parse_scanformatches(input_file: str):
    """
    Docstring pour parse_scanformatches
    Extrait les informations des matches à partir d'un fichier scan_for_matches.
    Args:
        input_file (str): Le chemin vers le fichier d'entrée scan_for_matches.
    Returns:
        list: Une liste de listes, chaque sous-liste contenant les champs d'un match.
    """
    seq_id = None
    positions = {}
    i=0
    with open(input_file, "r") as f:
        for line_traitee in f:
            line = line_traitee.strip()  # enlever les espaces en début/fin de ligne

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
    Docstring pour write_gff3
    Ecrire un fichier GFF3 à partir de la sortie de scan_for_matches
    :param input_file: fichier de sortie scan_for_matches
    :param output_file: fichier de sortie GFF3
    :param seq_fasta: fichier fasta optionnel pour obtenir la taille de la séquence si disponible

    return: message de confirmation
    """
    seq_id , positions = parse_scanformatches(input_file)
    source = 'scan_for_matches'
    taille = None
    feature_type = feature_type  # ou RBS/promoter/terminator selon le contexte
    if seq_fasta: # optionnel 
        taille = taille_seq(seq_fasta)

    with open(output_file, 'w') as out_fh:
        out_fh.write("##gff-version 3\n")

        if taille: # optionnel si taille connue 
            out_fh.write(f'##sequence-region {seq_id} 1 {taille}\n') # on indique la taille de la séquence
        
        for i  in sorted(positions):
            codon_start = None
            if len(positions[i]) == 3:
                start , end, codon_start= positions[i]
            else:
                start , end= positions[i]


            gff_start = min(int(start), int(end))
            gff_end = max(int(start), int(end)) # start <= end obligatoire en GFF3 d'apres la documentation

            score = '.' # pas de score dans scan_for_matches
            strand = '+' if int(start) < int(end) else '-'

            phase = '.' # pas de phase pour RBS/promoteur/terminateur
            if codon_start:
                attributs = f"ID{feature_type}={i};Note={source}_prediction;Start_Codon={codon_start}"
            else:
                attributs = f"ID{feature_type}={i};Note={source}_prediction"
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

