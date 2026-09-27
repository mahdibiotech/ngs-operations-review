# NGS Lot Review — prototype de revue opérationnelle

Projet personnel de **Mahdi Attabi** pour illustrer la rigueur d'exploitation bioinformatique. Données, taxons, seuils et référence entièrement **synthétiques**. Sans accès aux données, procédures, logiciels ou critères de PathoQuest. Non destiné à une décision clinique, BPF ou à une libération de lot.

## Le problème exploré

Après qu'une analyse NGS a produit une table de correspondances, comment faire passer **tout un lot** par les mêmes vérifications avant la revue humaine ? Un signal peut être accompagné de similarité à l'hôte, d'un soutien limité à une région ou d'un faible bruit dans le témoin négatif. L'outil conserve ces indices au lieu de conclure à une contamination ou d'éliminer automatiquement le signal.

PathoQuest décrit publiquement des essais de sécurité virale et des analyses donnant lieu à des rapports, dans un environnement de qualité réglementé. Ce projet explore une **interface de revue** déduite de ce contexte ; il ne présume d'aucune difficulté interne ni ne reproduit la plateforme iDTECT®. Voir `docs/requirements.md` pour les sources et les exigences.

## Lancer la démonstration (Python 3, bibliothèque standard)

```bash
python3 screen.py
python3 validate_demo.py
python3 -m unittest discover -s . -v
```

La première commande écrit :

- `results/report.json` : état du lot, résultats complets, motifs de revue et empreintes des entrées/règles ;
- `results/report.html` : rapport autonome à ouvrir dans un navigateur ;
- `results/worklist.tsv` : candidats à examiner par une personne.
- `results/validation_report.json` (deuxième commande) : huit défis synthétiques, résultat observé et empreintes des entrées.

Chaque correspondance est confrontée au catalogue versionné `config/synthetic_references.json` : accession, taxon et longueur doivent correspondre. Le nom de version doit correspondre à `metadata.reference_snapshot`. L'empreinte SHA-256 du catalogue figure dans le rapport : une modification des références est visible lors de la revue. `validate_demo.py` vérifie les témoins, le blocage, une référence inconnue, la répétabilité et le comportement exactement au seuil et juste en dessous. Une erreur donne un code retour non nul.

Dans WSL : `explorer.exe results` puis ouvrir `report.html`.

**Résultat attendu du lot nominal** : `REVIEW_READY`, deux candidats. SyntheticVirus-A dans TEST001 passe les seuils fictifs sans indicateur complémentaire. SyntheticVirus-C dans TEST002 passe les seuils, mais porte les motifs `NEGATIVE_BACKGROUND`, `LIMITED_REGION_SUPPORT` et `HOST_SIMILARITY` ; il reste en revue humaine.

## Tester les échecs des témoins

```bash
python3 screen.py --hits data/scenarios/negative_failure/hits.tsv --output /tmp/pq-negative.json --html /tmp/pq-negative.html --worklist /tmp/pq-negative.tsv
python3 screen.py --hits data/scenarios/positive_failure/hits.tsv --output /tmp/pq-positive.json --html /tmp/pq-positive.html --worklist /tmp/pq-positive.tsv
python3 screen.py --metadata data/scenarios/low_depth/metadata.json --output /tmp/pq-depth.json --html /tmp/pq-depth.html --worklist /tmp/pq-depth.tsv
```

Ces trois commandes retournent le code `2` avec `QC_BLOCKED` et une liste de travail en état `ON_HOLD_QC`. Ne pas enchaîner avec `&&` si vous souhaitez voir les trois cas. `docs/requirements.md` lie chaque exigence à un test.

## Périmètre exact

Le prototype lit des **preuves résumées** venant d'une étape amont fictive. Il n'analyse pas de FASTQ, ne réalise ni alignement, ni BLAST, ni assemblage, ni identification virale, ni estimation de sensibilité/spécificité. Les contrôles et seuils sont pédagogiques ; le HTML n'est pas un certificat d'analyse. Une transposition en contexte BPF exigerait une validation et une infrastructure qualité distinctes.

Ce banc d'essais prouve seulement que **le code de démonstration** réagit comme prévu à ces huit jeux fictifs ; il ne constitue ni une validation de méthode ni une qualification d'environnement réglementé. Pour exécuter sur d'autres entrées résumées, utiliser `--hits`, `--metadata`, `--config`, `--references` avec `screen.py`.
