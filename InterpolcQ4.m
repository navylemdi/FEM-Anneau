function N=InterpolcQ4(ksi)
% Cette fonction évalue numériquement la valeur des fonctions N aux
% coordonnées demandées en Ksi:
c=0.5;
N =(c)*[ (1-ksi),(1+ksi) ];