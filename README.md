# NGS Lot Review — version 5

Projet personnel de Mahdi Attabi : démonstrateur reproductible d'un petit flux de dépistage viral. Il ajoute à la v4 un **panel FASTA public**, un alignement Bowtie 2, une recherche BLAST+ des lectures et contigs, un assemblage SPAdes optionnel, une identification limitée aux références incluses, et des estimations de performance comparées à des étiquettes de vérité terrain.

> **Périmètre :** logiciel de démonstration sur un panel minuscule, et challenge reads simulés depuis les références incluses. Ce n'est ni un certificat d'analyse, ni un rapport clinique, ni un résultat de test sur patient/produit, ni une décision de libération, ni une méthode validée/qualifiée en BPF. Les valeurs de sensibilité/spécificité ne décrivent que le panel simulé exécuté et ne mesurent pas les performances cliniques ou industrielles.

## Références publiques

Le fichier `data/public_references/panel.fasta` contient deux références virales RefSeq et un leurre mitochondrial humain RefSeq. Accession, rôle, organisme, longueur, URL NCBI et SHA-256 figurent dans `config/public_reference_metadata.json`. La source reste identifiée dans les en-têtes FASTA. Le panel est volontairement petit et ne peut pas couvrir la diversité virale ou permettre d'interpréter une absence de hit comme une absence de virus.

Pour actualiser les fichiers depuis NCBI par accession, lancer `python3 fetch_public_references.py`, puis vérifier les nouveaux SHA-256 consignés dans le fichier de métadonnées avant une analyse.

| Accession | Séquence | Rôle |
|---|---|---|
| NC_045512.2 | SARS-CoV-2 Wuhan-Hu-1 | Référence virale |
| NC_001803.1 | HRSV A isolate UK/S2.ts1C/1995 | Référence virale |
| NC_012920.1 | Mitochondrie humaine | Leurre de séquence hôte |

## Installer les outils

Avec Conda ou Mamba :

```bash
conda env create -f environment.yml
conda activate ngs-lot-review-v5
bowtie2 --version
blastn -version
spades.py --version
```

## Créer un challenge panel étiqueté

Depuis la racine du dépôt, générer une série déterministe de lectures simulées, positives à plusieurs proportions virales et négatives à séquences mitochondriales :

```bash
python3 simulate_challenge_v5.py --outdir data/v5_challenge
```

Cela crée 27 échantillons par défaut : 24 positifs simulés (2 virus × 4 fractions × 3 réplicats) et 3 négatifs simulés. Les lectures et vérités terrain sont générées à partir des références locales ; elles ne sont pas des données publiques d'échantillons biologiques. Les options de taille, de fraction, de longueur et de graine sont visibles avec `--help`.

## Exécuter l'alignement et le BLAST

```bash
python3 pipeline_v5.py \
  --manifest data/v5_challenge/manifest.csv \
  --outdir results/v5
```

Le programme construit les index locaux, aligne les reads avec Bowtie 2, cherche les reads contre le FASTA avec BLAST+, mesure les reads alignés et la couverture de référence, puis demande un **soutien indépendant** de BLAST pour proposer un candidat. Il écrit `results/v5/report.json` et `results/v5/report.html`. Les commandes, logs, fichiers SAM et intermédiaires BLAST sont conservés sous `results/v5/work/` et `results/v5/logs/` pour permettre l'inspection.

## Ajouter l'assemblage

```bash
python3 pipeline_v5.py \
  --manifest data/v5_challenge/manifest.csv \
  --outdir results/v5_assembly \
  --assemble --assemble-sample POS_NC_045512_2_F0.2_R1
```

Avec `--assemble`, SPAdes en mode `--rnaviral` assemble les échantillons sélectionnés ; si aucun `--assemble-sample` n'est fourni, il assemble tous les échantillons. Les contigs sont ensuite interrogés par BLAST+ sur le même panel. Les contigs sont une preuve d'appui à examiner, pas une confirmation autonome d'identité.

## Sensibilité et spécificité

Le rapport compare les étiquettes du manifeste au résultat binaire du dépistage, présente TP/FN/TN/FP, les effectifs et les intervalles de confiance de Wilson à 95 %. Pour les positifs, la sensibilité exige la détection du taxon attendu ; pour les négatifs, une détection de n'importe quel taxon viral est un faux positif. Les échantillons `unknown` sont exclus. Des effectifs réduits produisent des intervalles larges. Une série de lectures simulées mesure le comportement du logiciel face à cette simulation ; elle ne permet pas d'estimer une sensibilité clinique, une limite de détection validée ou la spécificité d'une méthode réglementée.

Les seuils de `config/pipeline_v5.json` sont des valeurs d'exploration. Les changer modifie les métriques du panel. Ils ne sont pas des critères PathoQuest, des seuils diagnostiques, des limites BPF ou une validation de méthode.

## Tests logiciels

```bash
python3 -m unittest discover -s . -v
```

Les tests couvrent les lecteurs FASTQ, le CIGAR SAM, les filtres BLAST et les calculs de métriques. Ils ne remplacent pas une validation des outils externes ni une validation de méthode.

## Nature du rapport HTML

Le HTML est un **rapport exploratoire généré automatiquement**. Il inclut un avertissement visible, la provenance et l'empreinte du panel, les résultats par échantillon, les métriques et leurs intervalles. Il ne comporte ni signature, ni revue indépendante, ni contrôle des versions approuvé, ni procédure qualité, ni audit trail immuable. Il ne doit pas être nommé ou remis comme certificat d'analyse.

## Références de documentation

- NCBI RefSeq `NC_045512.2`: https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2
- NCBI RefSeq `NC_001803.1`: https://www.ncbi.nlm.nih.gov/nuccore/NC_001803.1
- NCBI RefSeq `NC_012920.1`: https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1
- NCBI BLAST+ User Manual: https://www.ncbi.nlm.nih.gov/books/NBK279690/
- Bowtie 2 Manual: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
- SPAdes documentation: https://ablab.github.io/spades/
- Pour la trajectoire d'évolution vers un environnement BPF et les preuves à produire, voir `docs/v5-method-and-limits.md`.
