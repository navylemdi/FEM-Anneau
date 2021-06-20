clear variables
close all

%% Calcul de convergence pour le T6
type_maille = 'T6';

% Listes des dimensions à tester
L_nRadius = [1,2,3,4,5,6,7,8,9];
L_nTheta = [1,2,3,4,5,6,7,8,9];

ind = 0;
indR=0;

for nRadius = L_nRadius
    indR = indR+1;
    indT=0;
    for nTheta = L_nTheta
        indT=indT+1;
        [uv,uv_th, tab_sigma, tab_sigma_th, Coord, Connect] = modele_T6(nRadius,nTheta);
        err_uv = erreur_depla(uv,uv_th,Coord,Connect,type_maille);
        err_tab_s = erreur_sigma(tab_sigma,tab_sigma_th,Coord,Connect,type_maille);
%         [indR,indT]
        E_uv(indR,indT) = err_uv;
        E_sx(indR,indT) = err_tab_s(1);
        E_sy(indR,indT) = err_tab_s(2);
        E_sxy(indR,indT) = err_tab_s(3);
        E_sr(indR,indT) = err_tab_s(4);
        E_st(indR,indT) = err_tab_s(5);
        E_srt(indR,indT) = err_tab_s(6);
    end
end




