function tab_err = erreur_sigma(tab_sigma,tab_sigma_th,Coord,Connect,type_maille)

    nelin=0;
    for i=1:size(Connect,1)
        if (Connect(i,4)~=0)
            nelin=nelin+1;
        end
    end

    tab_err = [0,0,0,0,0,0];
    if strcmp(type_maille,'T6')
        nnel = 6;
        % Paramètres pour l'intégration de Gauss
        npg = 6;
        [xloc,W]= GaussTriangle(npg);%Ou Quadrangle
        for ie=1:nelin
            nd = zeros(1,nnel);
            coor = zeros(nnel,2);
            % Extraire les infos propres à cet élément (coord., connect.,propriétés. ...)
            % Boucle sur les noeuds de l'élément en question
            for i=1:nnel
                nd(i) = Connect(ie,i);        % Connectivité des noeuds de l'élément
                coor(i,:) = Coord(nd(i),:);  % Coordonnées des noeuds de l'élément
                tab_sigma_elt(i,:) = tab_sigma(nd(i),:);
                tab_sigma_th_elt(i,:) = tab_sigma_th(nd(i),:);
            end
            % --------- Boucle sur les points d'intégration ----------
            for ig=1:npg
                % Extraction des coordonnees des points d'intégration
                k = xloc(1,ig);  % Ksi du point d'intégration sur le carré en cours
                n = xloc(2,ig);  % Eta du point d'intégration sur le carré en cours
                % Évaluation numérique de la matrice Bksieta pour le point et l'élément
                Bke= BKsiEtaT6(k,n);
                % Évaluation numérique des fonctions N pour le point et l'élément
                Nnum = InterpolT6(k,n);
                % Calcul numérique de la matrice Jacobienne (élément isoparamétrique)
                J= Bke*coor;
                % Calcul de son déterminant
                detJ= det(J);
                for i = 1:6
                    tab_err(i) = tab_err(i) + W(ig)*(Nnum*(tab_sigma_elt(:,i)-tab_sigma_th_elt(:,i)))^2*detJ;
                end
            end
        end
    end
    if strcmp(type_maille,'Q8')
        nnel = 8;
        % Paramètres pour l'intégration de Gauss
        npg = 9;
        [xloc,W]= GaussQuadrangle(npg);%Ou Quadrangle
        for ie=1:nelin
            nd = zeros(1,nnel);
            coor = zeros(nnel,2);
            % Extraire les infos propres à cet élément (coord., connect.,propriétés. ...)
            % Boucle sur les noeuds de l'élément en question
            for i=1:nnel
                nd(i) = Connect(ie,i);        % Connectivité des noeuds de l'élément
                coor(i,:) = Coord(nd(i),:);  % Coordonnées des noeuds de l'élément
                tab_sigma_elt(i,:) = tab_sigma(nd(i),:);
                tab_sigma_th_elt(i,:) = tab_sigma_th(nd(i),:);
            end
            % --------- Boucle sur les points d'intégration ----------
            for ig=1:npg
                % Extraction des coordonnees des points d'intégration
                k = xloc(ig,1);  % Ksi du point d'intégration sur le carré en cours
                n = xloc(ig,2);  % Eta du point d'intégration sur le carré en cours
                % Évaluation numérique de la matrice Bksieta pour le point et l'élément
                Bke= BKsiEtaQ8(k,n);
                % Évaluation numérique des fonctions N pour le point et l'élément
                Nnum = InterpolQ8(k,n);
                % Calcul numérique de la matrice Jacobienne (élément isoparamétrique)
                J= Bke*coor;
                % Calcul de son déterminant
                detJ= det(J);
                for i = 1:6
                    tab_err(i) = tab_err(i) + W(ig)*(Nnum*(tab_sigma_elt(:,i)-tab_sigma_th_elt(:,i)))^2*detJ;
                end
            end
        end
    end
    if strcmp(type_maille,'Q9')
        nnel = 9;
        
        % Paramètres pour l'intégration de Gauss
        npg = 9;
        [xloc,W]= GaussQuadrangle(npg);
        for ie=1:nelin
            nd = zeros(1,nnel);
            coor = zeros(nnel,2);
            % Extraire les infos propres à cet élément (coord., connect.,propriétés. ...)
            % Boucle sur les noeuds de l'élément en question
            for i=1:nnel
                nd(i) = Connect(ie,i);        % Connectivité des noeuds de l'élément
                coor(i,:) = Coord(nd(i),:);  % Coordonnées des noeuds de l'élément
                tab_sigma_elt(i,:) = tab_sigma(nd(i),:);
                tab_sigma_th_elt(i,:) = tab_sigma_th(nd(i),:);
            end
            % --------- Boucle sur les points d'intégration ----------
            for ig=1:npg
                % Extraction des coordonnees des points d'intégration
                k = xloc(ig,1);  % Ksi du point d'intégration sur le carré en cours
                n = xloc(ig,2);  % Eta du point d'intégration sur le carré en cours
                % Évaluation numérique de la matrice Bksieta pour le point et l'élément
                Bke= BKsiEtaQ9(k,n);
                % Évaluation numérique des fonctions N pour le point et l'élément
                Nnum = InterpolQ9(k,n);
                % Calcul numérique de la matrice Jacobienne (élément isoparamétrique)
                J= Bke*coor;
                % Calcul de son déterminant
                detJ= det(J);
                for i = 1:6
                    tab_err(i) = tab_err(i) + W(ig)*(Nnum*(tab_sigma_elt(:,i)-tab_sigma_th_elt(:,i)))^2*detJ;
                end
            end
        end
    end
    tab_err = sqrt(tab_err);
end
        