MathProj INSA 2024

Sujet: Optimisation numérique. Groupe I
Ylann Malherbe, Faiza Manzoor, Octavio Alonso Massa Perez

Le projet est constitué de deux fichiers différents: main.py et calcutil.py

main.py contient le code spécifique à chaque méthode et calcutil contient des fonctions de calculs utilitaires qu'on utilise plusieures fois.
On a codé les méthodes suivantes:

    - Méthode gradient à pas fixe (fixedStepGradientMethod)
    - Méthode gradient à pas optimal, avec méthode Wolfe pour trouver le bon pas (optimalStepGradientMethod)
    - Méthode Newton (newtonMethod)
    - Méthode BFGS-Newton (quasiNewtonMethod)

En exécutant main.py, le programme va selectionner des valeurs de r, h et lambda aléatoirement pour après converger vers le minimum on utilisant les différentes méthodes, pour après afficher les résultats de chaque méthode.


Dépendances:

    - Random (Bibliothèque standard)
    - Math (Bibliiothèque standard)


Exécution:

Pour exécuter le code il faut avoir installé Python 3.12.x et utiliser la commande "python main.py" dans le directeur racine du projet.

Pour modifier le volume V0, il faut juste modifier sa valeur dans la ligne 4 du fichier calcutil.py.

Pour modifier les valeurs de r, h et lambda, il faut juste modifier ces valeurs dans la fonction main du à la fin du fichier main.py