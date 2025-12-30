#!/usr/bin/env python3
# author: Laroussi Bachri
# December 29th , 2025
# M1 BBS , Université de Toulouse
# Projet Bioinformatique pour la génomique
import re
import sys

"""
Module pour parser la sortie de GeneMark et générer un fichier GFF3.
Le format GFF3 est décrit ici :
https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
Nous avons utilisé cette ressource pour construire le fichier GFF3.
"""

def extract_cds_genemark(input_file : str):
    """
    Docstring pour extract_cds_genemark
    Extrait les informations des CDS à partir d'un fichier GeneMark.
    Args:
        input_file (str): Le chemin vers le fichier d'entrée GeneMark.
    Returns:
        list: Une liste de listes, chaque sous-liste contenant les champs d'un CDS.
    """

    cds_list = []
    state = "outside" # état initial : en dehors de la section CDS , dans quelle section on est

    with open(input_file, "r") as f:
        for line_traitee in f:
            line = line_traitee.strip()  # enlever les espaces en début/fin de ligne

            if state == "outside":  # on cherche l'entête des CDS
                if "List of Open reading frames" in line:
                    state = "header" # trouvé 
                continue

            if state == "header":
                # on attend la "barre" sous l'entête de colonnes
                if line.startswith("--------"):
                    state = "data"
                continue

            # state == "data"
            if "List of Regions of interest" in line or "ABOUT THE MATRIX USED" in line: # fin des CDS
                break # sortir de la boucle pas d'informations utiles après

            if line == "" or line.startswith("Left") or line.startswith("end") or line.startswith("--------"): # ignorer lignes vides ou sans intérêt
                continue 

            # Les lignes CDS commencent par un nombre (après strip) et contiennent direct/complement + fr
            if line[0].isdigit() and (("direct" in line) or ("complement" in line)) and ("fr" in line):
                fields = line.split()  # on split par défaut sur les espaces , obtention sous forme de liste
                # fields attendu: [left, right, strand, 'fr', frame, coding_prob, start_prob]
                # Ex: ['2402','3313','direct','fr','2','0.47','0.06']
                cds_list.append(fields)
            

    print(cds_list)
    # choisir un CDS par (strand, frame, stop)

    groupes = {}

    for fields in cds_list:
        frame = int(fields[4])  # fr 1 -> 0, fr 2 -> 1, fr3 -> 2

        if fields[2] == "direct": # brin direct le stop est right end
            strand = "+"
            stop = int(fields[1])   # right end
        elif fields[2] == "complement": # brin complémentaire le stop est left end
            strand = "-"
            stop = int(fields[0])   # left end

        cle = (strand, frame, stop) # clé unique pour chaque 'groupe' partageant ces caractéristiques

        if cle not in groupes: # pas encore dans les groupes, on ajoute
            groupes[cle] = fields 
    return list(groupes.values()) # retourner la liste des CDS uniques

def parse_gm_file(input_file : str):
    """
    Docstring pour parse_gm_file
    Extrait l'identifiant et nom du programme utilisé.

    Args:
        input_file (str): Le chemin vers le fichier d'entrée GeneMark.
    return:

    seqid :  ID de la séquence (ex : 1_5404 ou nom du fragment FASTA)

    source :  Nom du programme (GeneMark, GeneMark.hmm, scan_for_matches)
    """
    seqid_pat = r'Sequence:\s(\S*)' # motif pour extraire l'ID de la séquence
    taille_pat = r'Sequence length:\s(\d+)' # motif pour extraire la longueur de la séquence
    id = ''
    source = 'Unknown' 
    taille = '.'
    with open(input_file , 'r') as fh:
        for line in fh:
            if line.startswith('Sequence:'):
                seqid_match = re.search(seqid_pat , line) 
                if seqid_match:
                    id = seqid_match.group(1) # extraire l'ID
                else:
                    id = 'unknown_id'
            if 'GENEMARK' in line.upper(): # détecter le programme GeneMark 
                # Le fichier ne presente pas explicitement le nom du programme utilisé
                source = 'GeneMark' # on suppose GeneMark si on trouve ce mot
            if line.startswith('Sequence length:'):
                taille_match = re.search(taille_pat , line)
                if taille_match:
                    taille = int(taille_match.group(1))
            if id and source != 'Unknown' and taille != '.':
                break # on a trouvé les deux informations nécessaires
    return id, source , taille
            
def write_gff3(input_file : str, output_file : str):
    """
    Docstring pour write_gff3
    Cette fonction écrit un fichier GFF3 à partir des données extraites d'un fichier GeneMark.

    Elle utilise les fonctions `parse_gm_file` et `extract_cds_genemark` pour obtenir les informations nécessaires.
    
    :param input_file: chemin du fichier GeneMark en entrée
    :param output_file: chemin du fichier GFF3 en sortie
    :return: message de confirmation de l'écriture du .GFF3
    """

    with open(output_file, 'w') as out_fh:
        # Ecrire l'en-tête GFF3
        seq_id , source , taille = parse_gm_file(input_file)
        cds_list = extract_cds_genemark(input_file)
        out_fh.write("##gff-version 3\n")
        out_fh.write(f'##sequence-region {seq_id} 1 {taille}\n')
        # Ecrire chaque CDS au format GFF3
        for idx, fields in enumerate(cds_list):
            start = fields[0] 
            end = fields[1] 
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