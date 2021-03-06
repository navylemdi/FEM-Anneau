function N=InterpolQ9(ksi,eta)
% Cette fonction évalue numériquement la valeur des fonctions N aux
% coordonnées demandées en Ksi et Eta:
N =[ (1-ksi)*(1-eta)*ksi*eta/4, ...        
     -(1-ksi^2)*(1-eta)*eta/2,  ...                    
    -(1+ksi)*(1-eta)*ksi*eta/4,    ...                                
    (1+ksi)*(1-eta^2)*ksi/2, ...
    (1+ksi)*(1+eta)*ksi*eta/4, ...
    (1-ksi^2)*(1+eta)*eta/2,...
    -(1-ksi)*(1+eta)*ksi*eta/4,...
    -(1-ksi)*(1-eta^2)*ksi/2,...
    (1-ksi^2)*(1-eta^2)];
end