# Prototype de revue opérationnelle — exigences et cahier d'essais

**Statut** : démonstrateur indépendant sur données entièrement synthétiques, non qualifié et non validé BPF. Les règles de `config/demo_rules.json` ne sont pas celles de PathoQuest.

## Besoin déduit des informations publiques

PathoQuest décrit une chaîne allant de la préparation et l'analyse bioinformatique à des rapports utilisables, avec des méthodes communes sur deux sites. Son offre Opérations demande exécution dans les délais, documentation, qualité et participation à la validation des outils. Ce démonstrateur explore un **point de passage de revue** après la génération de résultats d'alignement, sans supposer l'existence d'un problème interne ni reproduire iDTECT®.

Sources :
- https://pathoquest.com/quality-regulatory/
- https://pathoquest.com/services/adventitious-virus-testing/
- https://www.sfbi.fr/emplois/offre/202608260716-cdi-ingenieur-bio-informatique-operations
- https://www.sfbi.fr/emplois/offre/202602020427-stage-stage-m1m2-bioinformatique-optimisation-et-validation-de-pipelines-ngs

## Contrat d'entrée

- `metadata.json` : identifiant du lot, contexte, référence figée et échantillons (`test`, `negative`, `positive`). Chaque positif déclare le taxon attendu.
- `hits.tsv` : résultats *déjà calculés* par une étape amont ; taxon, référence, reads, régions distinctes, bases couvertes, taille de référence, identité et similarité à une séquence hôte. Le champ de similarité est fourni par l'amont dans cette fiction ; le programme ne réalise ni alignement ni recherche contre un génome hôte.
- `demo_rules.json` : valeurs illustratives, versionnées avec le programme. Aucun de ces nombres n'est un seuil analytique réel.

## Exigences et preuves

| ID | Attendu | Scénario de vérification |
| --- | --- | --- |
| R01 | Refuser les métriques impossibles et doublons échantillon/référence | `test_duplicate_reference_and_nonfinite_metrics_rejected` |
| R02 | Vérifier que **chaque** positif détecte son taxon attendu | `test_unexpected_positive_taxon_cannot_satisfy_control` |
| R03 | Mettre le lot en attente si le négatif atteint le seuil illustratif | `test_negative_control_failure_blocks_all_candidate_review` |
| R04 | Mettre le lot en attente si un échantillon est sous la profondeur minimale | `test_insufficient_depth_blocks_lot` |
| R05 | Garder un candidat ambigu visible et indiquer ses motifs | `test_normal_run_preserves_ambiguous_signal_for_review` |
| R06 | Lier les résultats aux entrées et règles exactes par SHA-256 | `test_provenance_changes_when_evidence_changes` |
| R07 | Échapper le texte externe dans le rapport HTML | `test_html_escapes_upstream_taxon` |

Codes de sortie : `0` = rapport produit, témoins du scénario passés ; `2` = rapport produit, lot bloqué ; `1` = entrée invalide ou erreur de lecture. **Aucun code ne signifie absence de virus, conformité d'un lot ou décision de libération.**

## Étapes nécessaires pour un vrai environnement BPF

Cadrer les formats et critères avec les équipes scientifique/qualité, jeux d'essais représentatifs, performances analytiques, versions de références, qualification informatique, revue et signatures, contrôles d'accès, sauvegarde, audit trail réglementaire, SOP approuvées et gestion des changements. Les empreintes SHA-256 du prototype ne constituent pas à elles seules une piste d'audit conforme.
# Extension v3 : références et défis reproductibles

## Extension FASTQ de démonstration

`fastq_demo.py` confronte quatre FASTQ fictifs à un FASTA fictif par correspondance exacte, sans tolérance aux erreurs. Les lectures ambiguës ne sont pas attribuées. La table produite est traitée par la même fonction `screen.run` que le scénario de preuves résumées. Le mode FASTQ conserve les empreintes des entrées et marque la similarité avec l'hôte comme non mesurée. Tests : `test_fastq_demo_reads_to_candidate_report` et `test_fastq_demo_rejects_malformed_quality_and_unknown_fasta`. Cette étape ne remplace pas une analyse NGS ou une validation de méthode.

| Exigence de démonstration | Preuve dans le dépôt | Limite |
| --- | --- | --- |
| Version de références cohérente avec le lot | `screen.py` contrôle `snapshot_id`, accession, taxon et longueur ; `test_reference_catalogue_rejects_mismatches_and_tracks_version` | Catalogue entièrement fictif, sans base de séquences. |
| Détecter une modification des références | `provenance.references_sha256` dans `results/report.json` | L'empreinte ne garantit pas l'intégrité d'une chaîne BPF. |
| Rejouer les comportements critiques | `python3 validate_demo.py` écrit `results/validation_report.json` ; `test_challenge_matrix_passes` | Huit cas synthétiques, pas une validation analytique. |

Ces contrôles illustrent les besoins de traçabilité et de répétabilité évoqués par les descriptions publiques de PathoQuest, sans supposer que l'entreprise rencontre ces problèmes en interne.
