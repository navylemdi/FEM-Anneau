function N=Interpolc(ksi)
%Cette fonction évalue numériquement la valeur des fonctions N aux
%coordonnées demandées en Ksi:
c=0.5;
N =(c)*[ -ksi*(1-ksi) ,2*(1-ksi^2) ,ksi*(1+ksi) ];
 