# V5 — méthode de démonstration et limites

## Ce que le pipeline exécute

1. Lit les FASTQ et le manifeste avec vérité terrain déclarée.
2. Construit un index Bowtie 2 depuis un panel RefSeq fixé par accession.
3. Aligne les lectures en local et calcule le nombre de reads ainsi que la fraction des positions virales couvertes.
4. Lance `blastn` contre une base nucléotidique locale construite depuis le même FASTA ; les identités et couvertures query sont filtrées par configuration.
5. Propose un candidat uniquement si le signal d'alignement, la couverture de la référence et le nombre de reads soutenus par BLAST dépassent les seuils de démonstration.
6. Sur demande, assemble les reads avec SPAdes `--rnaviral`, puis BLAST les contigs contre le panel.
7. Compare les candidats aux étiquettes connues, calcule sensibilité/spécificité au niveau échantillon et intervalles de Wilson, puis produit JSON et HTML.

## Ce que ces nombres signifient

Les fractions simulées sont des proportions de reads générées par le script, et non des concentrations virales mesurées dans un matériau. Le simulateur crée une petite diversité de substitution mais ne représente pas toute la variabilité des matrices, l'extraction, la préparation de librairie, les biais de séquençage, les contaminations croisées, les virus apparentés, les lots de réactifs ou les erreurs de vérité terrain.

Les métriques sont donc **des mesures de performance informatique sur un panel synthétique défini**. Elles ne sont pas une estimation transférable de la sensibilité/spécificité clinique ou industrielle. L'intervalle de Wilson aide à montrer l'incertitude d'échantillonnage ; il ne corrige pas les biais de conception du panel.

## Portée de l'identification

L'outil peut identifier uniquement les deux séquences virales fournies. Une absence d'alignement ou de BLAST peut venir d'un virus absent du panel, d'une divergence, d'une qualité faible ou de paramètres. Le leurre hôte est une séquence mitochondriale et ne représente pas le génome humain complet. Le rapport ne revendique donc ni dépistage virologique exhaustif ni exclusion d'agents viraux.

## Pourquoi ce n'est pas un certificat d'analyse

Un certificat d'analyse ou document réglementaire suppose un résultat traçable et autorisé selon un système qualité, avec méthode applicable, critères d'acceptation approuvés, vérification des données, revue habilitée, gestion des écarts, statut des instruments/logiciels, signatures et conservation contrôlée. Ce prototype ne fournit pas ces garanties. Le HTML est un rapport technique de démonstration, et doit être présenté comme tel.

## Éléments à construire avant toute revendication réglementaire

- Définir l'usage prévu, les matrices, la gamme d'agents, les critères d'acceptation et les responsabilités.
- Concevoir un protocole de validation avec matériaux de référence, positifs proches de la limite, négatifs représentatifs, virus apparentés et interférents, répétabilité/reproductibilité, robustesse, inclusivité/exclusivité, contamination croisée et critères statistiques préspécifiés.
- Tester des données et matériaux indépendants de ceux utilisés pour développer les paramètres ; documenter plan, résultats, déviations et approbation.
- Qualifier l'infrastructure, les versions, les dépendances, les droits, les sauvegardes, l'horodatage, la traçabilité, le contrôle des changements et la restauration.
- Ajouter authentification et signatures contrôlées, revue humaine, audit trail inviolable et génération approuvée des documents qualité.

Ces points constituent un plan de travail, pas la preuve que l'outil satisfait une norme ou qu'il est qualifié.
