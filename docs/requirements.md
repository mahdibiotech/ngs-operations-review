# Mini-spécification et cahier de tests — démonstration

**Objet.** Produire une revue reproductible de résultats d'alignement viral déjà calculés, sur données synthétiques. Le programme ne réalise ni mapping, ni assemblage, ni identification confirmatoire.

| ID | Exigence démontrée | Vérification |
| --- | --- | --- |
| F01 | Refuser les métriques invalides et les échantillons inconnus | Valeur négative ou identifiant inconnu → erreur |
| F02 | Signaler un hit selon les seuils illustratifs | 3 reads, couverture 5 %, identité 90 % → flag |
| F03 | Bloquer la revue si contrôle positif absent, contrôle négatif signalé ou profondeur insuffisante | Cas synthétiques en tests |
| F04 | Conserver empreinte SHA-256 des deux entrées et version du script | Comparer empreintes et fichiers |
| F05 | Mentionner clairement les limites d'interprétation | Champ limitations dans JSON |

**Hors périmètre.** Qualification logicielle, validation analytique, intégrité des données en environnement réglementé, gestion d'accès, LIMS, évaluation de la sensibilité/spécificité, taxonomie de référence, recherche d'homologie secondaire et décision de libération de lot. Une transposition BPF nécessiterait SOP approuvées, gestion des versions et changements, jeux de validation représentatifs, traçabilité encadrée et revue qualité.
