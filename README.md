# Acoustique_Matricant_Package
Ce package permet de realiser des calculs avec la méthode du matricant.
Le package cotient 5 classes :
- Proprietes_Physiques qui permet de récupérer les données physiques préalablement mises dans un fichier yaml;
- Matricant qui permet de réaliser les calculs jusqu'au matricant;
- Ploter qui permet d'afficher des figures des courbes ou des cmap en fonction des données d'entrée;
- RT_fluid qui calcul le coefficient de transmission pour un muliticouhe immergé;
- RT_solide qui permet de calculer le coefficient de transmission pour un multicouche compris entre 2 solides semi-infinis.

Un dossier exemple contient un fichier principale main.py et un un fichier yaml qui permettront de faire un calcul en incidence normale pour un multicouche immergé dans un fluide. Ces résultats sont comparés avec un cas de validation. Il y a un message d'erreur "invalid value encountered in arccos" dans la solution théorique, il est normal puisque physiquement il y a des bandes interdites.Je ne sais pas pourquoi ma fenêtre matplotlib ne s'affiche pas sur Pycharm, mais sur Spyder il n'y a pas de problème.

Le programme est vectorisé et utilise grandement le broadcasting. Ainsi, il faudrait prêter un grande attention au shape des élémements.
Le format des shape est le suivant :
- Nombre_couches    : nombre total de couche
- Nombre_pulsations : nombre total de pulsation
- Nombre_angles     : nombre d'angle
- dimension         : dimensions sptiales

Ce qui nous donne au final des shape du type (Nombre_couches,Nombre_pulsations,Nombre_angles,dimension,dimension). Vous comprendrez ainsi que même si vous utilisez un seul angle il faudra le définir comme un liste.

La méthode que j'ai utilisée est présentée dans la thèse de Cécile Baron "Le Développement en Série de Peano du Matricant Pour l’Etude de la Propagation des Ondes Elastiques en Milieux à Propriétés Continûment Variables" (2006).
