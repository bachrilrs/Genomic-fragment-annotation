# Projet d’annotation génomique

Annotation d’un fragment génomique bactérien
Auteur : Laroussi Bachri – M1 BBS, Université de Toulouse

## Description générale

Ce projet a pour objectif l’annotation structurale et fonctionnelle d’un fragment génomique bactérien (Bacillota), en utilisant les approches vues en cours de bioinformatique pour la génomique.

L’annotation structurale porte sur la prédiction des CDS, des signaux de régulation (promoteurs, RBS, terminateurs), et une annotation fonctionnelle qui est réalisée à partir des produits protéiques prédits.

Le language de programmation choisi est Python.
Le projet comprend également le développement de parsers en Python permettant de convertir les sorties de GeneMark, GeneMark.hmm et scan_for_matches vers le format standard GFF3.

⸻

## Organisation du projet

**Rapport_projet-Laroussi-Bachri.pdf**

Rapport final décrivant la démarche, les méthodes utilisées, l’analyse des résultats, les choix d’annotation et la synthèse biologique.

⸻

data/

Fichiers de sortie bruts des outils de prédiction :

- GeneMark (seuils 0.4 et 0.5)
- GeneMark.hmm
- Séquence génomique analysée

⸻

cds/

Séquences nucléotidiques des CDS retenues après annotation structurale, extraites à partir du fragment génomique.

⸻

prots/

Séquences protéiques correspondantes aux CDS, obtenues par traduction (EMBOSS transeq) et utilisées pour l’annotation fonctionnelle.

⸻

pat/

Motifs et matrices utilisés pour la recherche de signaux régulateurs avec scan_for_matches :

- RBS (strict, mismatch, matrices)
- Promoteurs (consensus, matrices)
- Terminateurs Rho-indépendants (tige-boucle)

⸻

results_sfm/

Résultats de scan_for_matches pour les différents motifs recherchés (RBS, promoteurs, terminateurs).

⸻

parser/

Scripts Python développés pour le projet :

- parser_gm.py : conversion GeneMark -> GFF3
- parser_genemarkhmm.py : conversion GeneMark.hmm -> GFF3
- parser_scanformatches.py : conversion scan_for_matches -> GFF3

Les pseudo-codes associés sont fournis dans des fichiers .txt séparés afin de préserver la lisibilité.

⸻

tables/

Tableaux de synthèse au format TSV :

- Comparaison des CDS entre outils
- Annotation structurale finale
- Annotation fonctionnelle
- Résultats intermédiaires (GeneMark, GeneMark.hmm, scan_for_matches)

⸻

Set-up

## Instructions pour reproduire l’annotation

### . Cloner le dépôt GitHub

```{bash}
git clone https://github.com/bachrilrs/Genomic-fragment-annotation
cd Genomic-fragment-annotation
tar -xvf scan_for_matches.tar
gcc -O -o scan_for_matches ggpunit.c scan_for_matches.c
```

Cette commande peut ne pas fonctionner sur certains systèmes. Le code source de scan_for_matches est ancien.

Essayer la commande suivante si la précédente échoue :

```{bash}
gcc -std=gnu89 -O -Wno-implicit-function-declaration -Wno-deprecated-non-prototype \
-o scan_for_matches ggpunit.c scan_for_matches.c -lm
```

### Tester scan_for_matches

```{bash}
chmod +x scan_for_matches
./run_tests tmp/
diff tmp test_output
```

### Executer les recherches de signaux régulateurs avec scan_for_matches en utilisant les fichiers dans le dossier `pat/`

```{bash}
./scan_for_matches pat/RBS_matrix_35 < data/seq1-5404.txt > results_sfm/RBS_matrix_35.out
````

## Remarques

- Les pseudo-codes sont fournis séparément du rapport principal pour garantir une meilleure lisibilité.
- Le format GFF3 est utilisé comme format d’échange standard pour l’annotation.
- Les choix d’annotation (bornes des CDS, chevauchements, signaux retenus) sont systématiquement justifiés dans le rapport.

⸻

## Outils principaux utilisés

- GeneMark / GeneMark.hmm
- scan_for_matches
- EMBOSS (extractseq, revseq, transeq)
- BLASTp
- SignalP, DeepTMHMM
- Python, R (pour la rédaction du projet)

## Compétences

- Python (parsing, CLI tools)
- Redaction d'un rapport scientifique
- Annotation génomique (CDS, RBS, promoters, terminators)
- Format de fichiers: FASTA, GFF3
- Pipelines
