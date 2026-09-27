# Extension sur données publiques : protocole de recherche

**Statut : protocole à exécuter, aucun résultat sur données réelles revendiqué.** L'environnement de cette démonstration ne dispose ni des exécutables bioinformatiques ni des FASTQ téléchargés. Ce document sépare les étapes programmables des preuves qu'il reste à produire.

## Exemple vérifiable

| Élément | Source publique | Usage et réserve |
| --- | --- | --- |
| Lectures Illumina `SRR13178806` | https://www.ncbi.nlm.nih.gov/sra/SRR13178806 | Cas positif attendu SARS-CoV-2 selon les métadonnées du dépôt ; environ 196 844 spots, 15,5 Mo de téléchargement indiqué par SRA. Ces métadonnées ne valent pas preuve indépendante pour une estimation de sensibilité. |
| Référence `MN908947.3` | https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3 | Référence SARS-CoV-2 ; conserver accession, version, fichier et SHA-256. |

Ne pas confondre cette séquence virale publique avec un panel d'épreuve de sécurité virale de biothérapies : la matrice, la préparation et l'usage diffèrent.

## Commandes de recherche sous Linux

Installer dans un environnement dédié des versions fixées et consignées de `sra-tools`, `minimap2`, `samtools`, `blast`, `spades`, `fastp` selon leurs notices officielles. Vérifier le type de lecture et le consentement/conditions de réutilisation du jeu de données. Disposer de l'espace disque nécessaire.

```bash
mkdir -p public_data/reference public_data/reads public_data/work
prefetch SRR13178806 --output-directory public_data/reads
fasterq-dump public_data/reads/SRR13178806/SRR13178806.sra \
  --split-files -O public_data/reads
efetch -db nucleotide -id MN908947.3 -format fasta > public_data/reference/MN908947.3.fasta
sha256sum public_data/reads/SRR13178806*.fastq public_data/reference/MN908947.3.fasta
```

Selon la structure de l'archive, `fasterq-dump` peut produire un ou deux FASTQ : inspecter les fichiers avant d'appliquer une commande paired-end. Conserver l'archive d'origine et les métadonnées SRA. Si un outil manque ou échoue, arrêter l'analyse et consigner l'erreur.

## Analyses à implémenter et vérifier sur ce jeu réel

1. **QC et préparation :** nombre de lectures, qualités par cycle, adaptateurs, lectures retenues et versions de `fastp`. Préétablir les critères de rejet et leur justification ; ne pas choisir les seuils après avoir vu le résultat.
2. **Alignement :** `minimap2 -ax sr` sur la référence fixée ; `samtools sort/index/flagstat/depth` pour produire nombre de lectures alignées, profondeur et couverture. Vérifier les lectures multimappées et les régions faiblement couvertes. Une cible virale connue ne suffit pas à mesurer la spécificité contre une base diversifiée.
3. **Assemblage :** assembler les lectures pertinentes avec `spades.py --only-assembler` et conserver les contigs, leur longueur, couverture et provenance. Pour un génome viral entier, adapter les paramètres à l'architecture de la bibliothèque et vérifier les artefacts d'amplicons.
4. **Identification :** créer une base locale `makeblastdb` à partir d'un panel versionné et lancer `blastn` sur les contigs ; conserver accession, identité, longueur d'alignement, couverture de requête, score et autres hits plausibles. Une correspondance BLAST isolée ne prouve ni l'identité biologique ni l'absence d'autres agents.
5. **Évaluation :** constituer un panel d'échantillons **indépendants** avec étiquettes positives/négatives vérifiées et matrice représentative ; préenregistrer un seuil et une définition au niveau échantillon ; calculer TP/FN/TN/FP, sensibilité = TP/(TP+FN), spécificité = TN/(TN+FP) et intervalles de confiance. Étudier séparément répétabilité, limite de détection, robustesse, contaminations croisées et domaines non couverts. Si TP+FN ou TN+FP vaut zéro, indiquer « non estimable ».

## Dossier BPF : éléments à faire approuver dans l'organisation

La validation analytique et la qualification du système demandent notamment un usage prévu défini, une analyse de risques, exigences approuvées, traçabilité exigences/tests, versions et changements contrôlés, données d'essai et critères d'acceptation préspécifiés, preuves d'exécution, gestion des anomalies, audit trail, accès et sauvegardes, qualification d'infrastructure et revue/approbation qualité. La portée dépend de l'usage et du système réellement installés. L'outil ne génère donc **pas** de certificat d'analyse : un tel document dépend d'un processus approuvé et d'une personne habilitée.

Références : [ICH Q2(R2)](https://database.ich.org/sites/default/files/ICH_Q2%28R2%29_Guideline_2023_1130.pdf), [EudraLex annexe 11](https://health.ec.europa.eu/system/files/2016-11/annex11_01-2011_en_0.pdf), [EudraLex annexe 15](https://health.ec.europa.eu/system/files/2016-11/2015-10_annex15_0.pdf), [NCBI BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/), [NCBI SRA Toolkit](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump).
