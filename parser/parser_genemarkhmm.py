#!/usr/bin/env python3
# author: Laroussi Bachri
# December 29th , 2025
# M1 BBS , Université de Toulouse
# Projet Bioinformatique pour la génomique
import re
import sys

"""
Module pour parser la sortie de GeneMark.hmm et générer un fichier GFF3.
Le format GFF3 est décrit ici :
https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
Nous avons utilisé cette ressource pour construire le fichier GFF3.
"""
def extract_infos_GenemarkHMM(input_file: str):
    """
    Docstring pour extract_infos_GenemarkHMM
    Extrait l'identifiant de la séquence, le nom du programme et les informations des
    CDS à partir d'un fichier de sortie genemarkhmm.
    :param input_file: chemin du fichier de sortie genemarkhmm
    :return: tuple (seq_id, source, cds_list)
        seq_id : ID de la séquence (ex : nom du fragment FASTA)
        source : Nom du programme (ex : GeneMark.hmm)
        cds_list : liste des CDS extraits du fichier
    """
    seq_id = ''
    source = ''
    cds_list = []
    with open(input_file , 'r') as fh:
        state = 'outside' # état initial avant la section des données
        for line_traitee in fh:
            line = line_traitee.strip()
            if 'Version' in line: # ligne contenant le nom du programme
                source = line.split()[0]
            if line.startswith('FASTA definition line:'): # ligne contenant l'ID de la séquence
                seq_match = re.search(r':\s(\S*)' , line)
                if seq_match:
                    seq_id = seq_match.group(1)
            if line.startswith('Predicted genes'): # début de la section des données
                state = 'data'
                continue

            if state == 'data' and line and line[0].isdigit(): # verifier que la ligne commence par un chiffre
                fields = line.split() # on recupère les champs
                cds_list.append(fields) 
            

    return seq_id , source , cds_list

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



def nettoyer_int(x: str):
    """
    Docstring pour nettoyer_int
    Eviter les erreurs de conversion d'une chaîne de caractères en entier en supprimant les >/<.
    (Notamment sur les positions dans genemarkhmm)
    :param x: string
    return: int
    """
    return int(x.replace("<","").replace(">",""))



def write_gff3(input_file: str , output_file: str ,seq_fasta=None):
    """
    Docstring pour write_gff3
    Ecrire un fichier GFF3 à partir de la sortie de genemarkhmm
    :param input_file: fichier de sortie genemarkhmm
    :param output_file: fichier de sortie GFF3
    :param seq_fasta: fichier fasta optionnel pour obtenir la taille de la séquence si disponible
    
    Cette fonction utilise les extract_infos_GenemarkHMM et taille_seq pour obtenir les informations nécessaires.

    return: message de confirmation
    """
    seq_id , source , cds_list= extract_infos_GenemarkHMM(input_file)
    taille = None
    if seq_fasta: # optionnel 
        taille = taille_seq(seq_fasta)

    with open(output_file, 'w') as out_fh:
        out_fh.write("##gff-version 3\n")

        if taille: # optionnel si taille connue 
            out_fh.write(f'##sequence-region {seq_id} 1 {taille}\n') # on indique la taille de la séquence
        
        for i , fields in enumerate(cds_list):
            left = nettoyer_int(fields[2]) 
            right = nettoyer_int(fields[3])
            start = min(left, right) # start <= end obligatoire en GFF3 d'apres la documentation
            end = max(left, right)

            score = '.' # pas de score dans genemarkhmm
            strand = fields[1] # + ou -

            if strand == '+': # brin direct
                phase = (start - 1) % 3 # on calcule la phase a partir de la position de début
            else:
                phase = (end - 1) % 3 # brin reverse on calcule a partir de la position de fin

            classe = fields[5] if len(fields)>=6 else '.' # entre typical et atypical (transfert horizontaux par ex)
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

    