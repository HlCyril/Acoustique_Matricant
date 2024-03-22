# Acoustique_Matricant_Package
Ce package permet de realiser des calculs avec la méthode du matricant.
Le package cotient 5 classes :
- Proprietes_Physiques qui permet de récupérer les données physiques préalablement mises dans un fichier yaml;
- Matricant qui permet de réaliser les calculs jusqu'au matricant;
- Ploter qui permet d'afficher des figures des courbes ou des cmap en fonction des données d'entrée;
- RT_fluid qui calcul le coefficient de transmission pour un muliticouhe immergé;
- RT_solide qui permet de calculer le coefficient de transmission pour un multicouche compris entre 2 solides semi-infinis.

Un dossier example contient un fichier principale main.py et un un fichier yaml qui permettront de faire un calcul en incidence normale pour un multicouche immergé dans un fluide. Ces résultats sont comparés avec un cas de validation.

Le programme est vectorisé et utilise grandement le broadcasting. Ainsi, il faudrait prêter un grande attention au shape des élémements.
Le format des shape est le suivant :
- Nombre_couches    : nombre total de couche
- Nombre_pulsations : nombre total de pulsation
- Nombre_angles     : nombre d'angle
- dimension         : dimensions sptiales
Ce qui nous donne au final des shape du type (Nombre_couches,Nombre_pulsations,Nombre_angles,dimension,dimension). Vous comprendrez ainsi que même si vous utilisez un seul angle il faudra le définir comme un liste.
