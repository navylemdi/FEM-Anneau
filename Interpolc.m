function N=Interpolc(ksi)
%Cette fonction �value num�riquement la valeur des fonctions N aux
%coordonn�es demand�es en Ksi:
c=0.5;
N =(c)*[ -ksi*(1-ksi) ,2*(1-ksi^2) ,ksi*(1+ksi) ];
 