function [GX,GY,GXY]=Gint(xloc,W,Nui,nnel,nelin,ddln,npg,D,matConnec,matNoeud,GX,GY,GXY,sol,dim)
if dim=='Q8'
    for ie=1:nelin  % Pour tous les éléments internes
        % Extraires les déplacements propres à cet élément à partir de ¨sol¨
        % Initialisation du vecteur colonne des déplacements pour cet élément
        uv=zeros(Nui,1);
        gex= zeros(nnel,1);
        gey= zeros(nnel,1);
        gexy= zeros(nnel,1);
        deform= zeros(3,1);
        contraintes= zeros(3,1);
        
        % Repérage de la position des noeuds de cet élément
        nd = matConnec(ie,:);  % Numéros des noeuds
        Kloc=Feeldof(nd,nnel,ddln); % Retourne les # de DDL concernés par l'élément
        % extraire le vecteur solution elementaire  des u et v pour les ¨nnel¨ noeuds de l'élément
        for l=1:Nui
            temp=Kloc(1,l);
            uv(l,1)=sol(temp,1);
        end
        % extraires les coordonnees
        for i=1:nnel
            nd(i) = matConnec(ie,i);        % Connectivité des noeuds de l'élément
            coor(i,:) = matNoeud(nd(i),:);  % Coordonnées des noeuds de l'élément
            
        end
        coor;
        % ------ Calcul des contraintes aux points d'intégration ----------
        % Boucle sur les points d'intégration
        for ig=1:npg
            % Extraction des coordonnees des points d'intégration
            k = xloc(ig,1); % Ksi du point d'intégration sur le carré en cours
            n = xloc(ig,2);  % Eta du point d'intégration sur le carré en cours
            % Évaluation numérique des fonctions N pour le point et l'élément
            Nnum = InterpolQ8(k,n);
            Bke= BKsiEtaQ8(k,n);
            xg=Nnum*coor(:,1);
            yg=Nnum*coor(:,2);
            % Calcul numérique de la matrice Jacobienne (élément isoparamétrique)
            J= Bke*coor;
            % Calcul de son inverse
            j= J^-1 ;
            % Calcul de son déterminant
            %(utile car l'intégration se fait sur élément de référence plus loin)
            DJ= det(J);
            Bxy = j*Bke;
            % Construction de la matrice B car élasticité 2D (avec les termes de Bxy)
            % RAPPEL: {epsilon}e=[B]*{U}e
            % [B](3*16) fait le lien entre les 16 déplacements et les 3 déformations
            % Initialisation
            B=zeros(3,Nui);
            % Assignation directe des 16 termes de B non nuls
            B(1,1)=Bxy(1,1);
            B(1,3)=Bxy(1,2);
            B(1,5)=Bxy(1,3);
            B(1,7)=Bxy(1,4);
            B(1,9)=Bxy(1,5);
            B(1,11)=Bxy(1,6);
            B(1,13)=Bxy(1,7);
            B(1,15)=Bxy(1,8);
            
            B(2,2)=Bxy(2,1);
            B(2,4)=Bxy(2,2);
            B(2,6)=Bxy(2,3);
            B(2,8)=Bxy(2,4);
            B(2,10)=Bxy(2,5);
            B(2,12)=Bxy(2,6);
            B(2,14)=Bxy(2,7);
            B(2,16)=Bxy(2,8);
            
            B(3,1)=Bxy(2,1);
            B(3,2)=Bxy(1,1);
            B(3,3)=Bxy(2,2);
            B(3,4)=Bxy(1,2);
            B(3,5)=Bxy(2,3);
            B(3,6)=Bxy(1,3);
            B(3,7)=Bxy(2,4);
            B(3,8)=Bxy(1,4);
            B(3,9)=Bxy(2,5);
            B(3,10)=Bxy(1,5);
            B(3,11)=Bxy(2,6);
            B(3,12)=Bxy(1,6);
            B(3,13)=Bxy(2,7);
            B(3,14)=Bxy(1,7);
            B(3,15)=Bxy(2,8);
            B(3,16)=Bxy(1,8);
            % Calcul et stockage des trois déformations pour le point en cours
            % RAPPEL: {epsilon}e=[B]*{U}e
            deform=B*uv;
            % Calcul et stockage des contraintes au point d'intégration en cours
            % RAPPEL: contraintes=[D]*{epsilon}e
            contraintes=D*deform;  % ICI ON a les contraintes de l'element ie et au point de Gauss IG
            gex= gex+W(ig)*(Nnum')*contraintes(1)*DJ;
            gey= gey+W(ig)*(Nnum')*contraintes(2)*DJ;
            gexy= gexy+W(ig)*(Nnum')*contraintes(3)*DJ;
        end    % fin sur les points d integration
        Kloc=Feeldof(nd,nnel,1);% Retourne les # de DDL concernés par l'élément
        GX= AssemblF(GX,gex,Kloc);% Assemble gex dans GX connaissant Kloc2
        GY= AssemblF(GY,gey,Kloc);% Assemble gey dans GY connaissant Kloc2
        GXY= AssemblF(GXY,gexy,Kloc);% Assemble gexy  dans GXY connaissant Kloc2
        
    end
elseif dim=='Q9'
    for ie=1:nelin  % Pour tous les éléments internes
        % Extraires les déplacements propres à cet élément à partir de ¨sol¨
        % Initialisation du vecteur colonne des déplacements pour cet élément
        uv=zeros(Nui,1);
        gex= zeros(nnel,1);
        gey= zeros(nnel,1);
        gexy= zeros(nnel,1);
        deform= zeros(3,1);
        contraintes= zeros(3,1);
        % Repérage de la position des noeuds de cet élément
        nd = matConnec(ie,:);  % Numéros des noeuds
        Kloc=Feeldof(nd,nnel,ddln); % Retourne les # de DDL concernés par l'élément
        % extraire le vecteur solution elementaire  des u et v pour les ¨nnel¨ noeuds de l'élément
        for l=1:Nui
            temp=Kloc(1,l);
            uv(l,1)=sol(temp,1);
        end
        % extraires les coordonnees
        for i=1:nnel
            nd(i) = matConnec(ie,i);        % Connectivité des noeuds de l'élément
            coor(i,:) = matNoeud(nd(i),:);  % Coordonnées des noeuds de l'élément
        end
        coor;
        nd;
        % ------ Calcul des contraintes aux points d'intégration ----------
        % Boucle sur les points d'intégration
        for ig=1:npg
            % Extraction des coordonnees des points d'intégration
            k = xloc(ig,1);  % Ksi du point d'intégration sur le carré en cours
            n = xloc(ig,2);  % Eta du point d'intégration sur le carré en cours
            % Évaluation numérique des fonctions N pour le point et l'élément
            Nnum = InterpolQ9(k,n);
            Bke= BKsiEtaQ9(k,n);
            %xg=Nnum*coor(:,2);
            %yg=Nnum*coor(:,1);
            % Calcul numérique de la matrice Jacobienne (élément isoparamétrique)
            J= Bke*coor;
            % Calcul de son inverse
            j= J^-1 ;
            % Calcul de son déterminant
            %(utile car l'intégration se fait sur élément de référence plus loin)
            DJ= det(J);
            Bxy = j*Bke;
            % Construction de la matrice B car élasticité 2D (avec les termes de Bxy)
            % RAPPEL: {epsilon}e=[B]*{U}e
            % [B](3*16) fait le lien entre les 16 déplacements et les 3 déformations
            % Initialisation
            B=zeros(3,nnel*2);
            % Assignation directe des 16 termes de B non nuls
            B(1,1)=Bxy(1,1);
            B(1,3)=Bxy(1,2);
            B(1,5)=Bxy(1,3);
            B(1,7)=Bxy(1,4);
            B(1,9)=Bxy(1,5);
            B(1,11)=Bxy(1,6);
            B(1,13)=Bxy(1,7);
            B(1,15)=Bxy(1,8);
            B(1,17)=Bxy(1,9);
            
            B(2,2)=Bxy(2,1);
            B(2,4)=Bxy(2,2);
            B(2,6)=Bxy(2,3);
            B(2,8)=Bxy(2,4);
            B(2,10)=Bxy(2,5);
            B(2,12)=Bxy(2,6);
            B(2,14)=Bxy(2,7);
            B(2,16)=Bxy(2,8);
            B(2,18)=Bxy(2,9);
            
            B(3,1)=Bxy(2,1);
            B(3,2)=Bxy(1,1);
            B(3,3)=Bxy(2,2);
            B(3,4)=Bxy(1,2);
            B(3,5)=Bxy(2,3);
            B(3,6)=Bxy(1,3);
            B(3,7)=Bxy(2,4);
            B(3,8)=Bxy(1,4);
            B(3,9)=Bxy(2,5);
            B(3,10)=Bxy(1,5);
            B(3,11)=Bxy(2,6);
            B(3,12)=Bxy(1,6);
            B(3,13)=Bxy(2,7);
            B(3,14)=Bxy(1,7);
            B(3,15)=Bxy(2,8);
            B(3,16)=Bxy(1,8);
            B(3,17)=Bxy(2,9);
            B(3,18)=Bxy(1,9);
            % Calcul et stockage des trois déformations pour le point en cours
            % RAPPEL: {epsilon}e=[B]*{U}e
            deform=B*uv;
            % Calcul et stockage des contraintes au point d'intégration en cours
            % RAPPEL: contraintes=[D]*{epsilon}e
            contraintes=D*deform;  % ICI ON a les contraintes de l'element ie et au point de Gauss IG
            gex= gex+W(ig)*(Nnum')*contraintes(1,1)*DJ;
            gey= gey+W(ig)*(Nnum')*contraintes(2,1)*DJ;
            gexy= gexy+W(ig)*(Nnum')*contraintes(3,1)*DJ;
        end    % fin sur les points d integration
        Kloc2=Feeldof(nd,nnel,1);% Retourne les # de DDL concernés par l'élément
        GX= AssemblF(GX,gex,Kloc2);% Assemble gex dans GX connaissant Kloc2
        GY= AssemblF(GY,gey,Kloc2);% Assemble gey dans GY connaissant Kloc2
        GXY= AssemblF(GXY,gexy,Kloc2);% Assemble gexy  dans GXY connaissant Kloc2
    end % Fin de la boucle sur les éléments
elseif dim=='T6'
    for ie=1:nelin  % Pour tous les éléments internes
        % Extraires les déplacements propres à cet élément à partir de ¨sol¨
        % Initialisation du vecteur colonne des déplacements pour cet élément
        uv=zeros(Nui,1);
        gex= zeros(nnel,1);
        gey= zeros(nnel,1);
        gexy= zeros(nnel,1);
        deform= zeros(3,1);
        contraintes= zeros(3,1);
        % Repérage de la position des noeuds de cet élément
        nd = matConnec(ie,:);  % Numéros des noeuds
        Kloc=Feeldof(nd,nnel,ddln); % Retourne les # de DDL concernés par l'élément
        % extraire le vecteur solution elementaire  des u et v pour les ¨nnel¨ noeuds de l'élément
        for l=1:Nui
            temp=Kloc(1,l);
            uv(l,1)=sol(temp,1);
        end
        % extraires les coordonnees
        for i=1:nnel
            nd(i) = matConnec(ie,i);        % Connectivité des noeuds de l'élément
            coor(i,:) = matNoeud(nd(i),:);  % Coordonnées des noeuds de l'élément
        end
        coor;
        nd;
        % ------ Calcul des contraintes aux points d'intégration ----------
        % Boucle sur les points d'intégration
        for ig=1:npg
            % Extraction des coordonnees des points d'intégration
            k = xloc(1,ig);  % Ksi du point d'intégration sur le carré en cours
            n = xloc(2,ig);  % Eta du point d'intégration sur le carré en cours
            % Évaluation numérique des fonctions N pour le point et l'élément
            Nnum = InterpolT6(k,n);
            Bke= BKsiEtaT6(k,n);
            %xg=Nnum*coor(:,1);
            %yg=Nnum*coor(:,2);
            %plot(xg,yg,'.b')
            % Calcul numérique de la matrice Jacobienne (élément isoparamétrique)
            J= Bke*coor;
            % Calcul de son inverse
            j= J^-1 ;
            % Calcul de son déterminant
            %(utile car l'intégration se fait sur élément de référence plus loin)
            DJ= det(J);
            Bxy = j*Bke;
            % Construction de la matrice B car élasticité 2D (avec les termes de Bxy)
            % RAPPEL: {epsilon}e=[B]*{U}e
            % [B](3*16) fait le lien entre les 16 déplacements et les 3 déformations
            % Initialisation
            B=zeros(3,Nui);
            % Assignation directe des 16 termes de B non nuls
            B(1,1)=Bxy(1,1);
            B(1,3)=Bxy(1,2);
            B(1,5)=Bxy(1,3);
            B(1,7)=Bxy(1,4);
            B(1,9)=Bxy(1,5);
            B(1,11)=Bxy(1,6);
            
            B(2,2)=Bxy(2,1);
            B(2,4)=Bxy(2,2);
            B(2,6)=Bxy(2,3);
            B(2,8)=Bxy(2,4);
            B(2,10)=Bxy(2,5);
            B(2,12)=Bxy(2,6);
            
            B(3,1)=Bxy(2,1);
            B(3,2)=Bxy(1,1);
            B(3,3)=Bxy(2,2);
            B(3,4)=Bxy(1,2);
            B(3,5)=Bxy(2,3);
            B(3,6)=Bxy(1,3);
            B(3,7)=Bxy(2,4);
            B(3,8)=Bxy(1,4);
            B(3,9)=Bxy(2,5);
            B(3,10)=Bxy(1,5);
            B(3,11)=Bxy(2,6);
            B(3,12)=Bxy(1,6);
            % Calcul et stockage des trois déformations pour le point en cours
            % RAPPEL: {epsilon}e=[B]*{U}e
            deform=B*uv;
            % Calcul et stockage des contraintes au point d'intégration en cours
            % RAPPEL: contraintes=[D]*{epsilon}e
            contraintes=D*deform;  % ICI ON a les contraintes de l'element ie et au point de Gauss IG
            gex= gex+W(ig)*(Nnum')*contraintes(1,1)*DJ;
            gey= gey+W(ig)*(Nnum')*contraintes(2,1)*DJ;
            gexy= gexy+W(ig)*(Nnum')*contraintes(3,1)*DJ;
            
        end    % fin sur les points d integration
        Kloc2=Feeldof(nd,nnel,1);% Retourne les # de DDL concernés par l'élément
        GX= AssemblF(GX,gex,Kloc2);% Assemble gex dans GX connaissant Kloc2
        GY= AssemblF(GY,gey,Kloc2);% Assemble gey dans GY connaissant Kloc2
        GXY= AssemblF(GXY,gexy,Kloc2);% Assemble gexy  dans GXY connaissant Kloc2
        
    end % Fin de la boucle sur les éléments
else
    disp('Element non pris en compte')
end
end