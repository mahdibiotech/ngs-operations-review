# NGS Operations Review — démonstration pour PathoQuest

Projet personnel de candidature de Mahdi Attabi. Données et noms de virus fictifs. Sans lien avec les données, logiciels ou procédures propriétaires de PathoQuest.

## Pourquoi ce projet
Ce projet associe exécution d'analyses NGS, contrôles, rendu et documentation. Cette petite démonstration illustre le traitement vérifiable d'un tableau de hits déjà produits : validation des entrées, revue des contrôles, signalements illustratifs, empreintes des fichiers et rapport JSON.

## Exécution (Python 3, bibliothèque standard)

```bash
python3 screen.py --hits data/hits.tsv --metadata data/metadata.json --output report.json
python3 -m unittest -v
```

Le code de sortie 2 signale un blocage QC, 1 une entrée invalide. Lire `docs/requirements.md` pour les exigences et les limites. Aucun seuil du script ne correspond à un critère PathoQuest ; aucun résultat ne permet de conclure à une contamination ou à une conformité BPF.

## Suite possible en poste
Après cadrage avec l'équipe : intégrer les formats réels validés et les contrôles requis, définir les critères avec les responsables scientifiques et qualité, documenter les cas d'échec, exécuter une validation formalisée, puis intégrer à l'environnement d'exploitation selon leurs procédures.
