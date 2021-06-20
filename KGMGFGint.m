function [KG,FG,MG]=KGMGFGint(xloc,W,Nui,nnel,nelin,ddln,npg,D,matConnec,matNoeud,KG,MG,FG,dim)
if dim=='Q8'
    for ie=1:nelin
        % Extraire les infos propres à cet élément (coord., connect.,propriétés. ...)
        % Boucle sur les noeuds de l'élément en question
        for i=1:nnel
            nd(i) = matConnec(ie,i);        % Connectivité des noeuds de l'élément
            coor(i,:) = matNoeud(nd(i),:);  % Coordonnées des noeuds de l'élément
        end
        coor;
        nd;
        % Initialisation de la matrice de rigidité élémentaire
        ke=zeros(Nui,Nui);
        me= zeros(nnel,nnel);
        % Initialisation du vecteur sollicitation élémentaire
        fe=zeros(Nui,1);
        % --------- Boucle sur les points d'intégration ----------
        for i=1:npg
            % Extraction des coordonnees des points d'intégration
            k = xloc(i,1);  % Ksi du point d'intégration sur le carré en cours
            n = xloc(i,2);  % Eta du point d'intégration sur le carré en cours
            % Évaluation numérique de la matrice Bksieta pour le point et l'élément
            Bke= BKsiEtaQ8(k,n);
            % Évaluation numérique des fonctions N pour le point et l'élément
            Nnum = InterpolQ8(k,n);
            xg=Nnum*coor(:,1);
            yg=Nnum*coor(:,2);
            % Calcul numérique de la matrice Jacobienne (élément isoparamétrique)
            J= Bke*coor;
            % Calcul de son inverse
            j= J^-1 ;
            % Calcul de son déterminant
            %(utile car l'intégration se fait sur élément de référence plus loin)
            DJ= det(J);
            % Calcul numérique de la matrice Bxy
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
            % Contruction numérique de la matrice de rigidite elementaire
            % Note: Sommation des (poids) X (fonction évaluée au point de Gauss)
            ke = ke  + W(i)*(B')*D*B*DJ;
            me = me  + W(i)*(Nnum')*Nnum*DJ;
            % Calcul de la matrice [N](2x16 dans le cas des Q8) à partir de Nnum
            % Initialisation
            N=zeros(2,nnel*2);
            % Assignation directe des 16 termes de N non nuls
            N(1,1)=Nnum(1);
            N(1,3)=Nnum(2);
            N(1,5)=Nnum(3);
            N(1,7)=Nnum(4);
            N(1,9)=Nnum(5);
            N(1,11)=Nnum(6);
            N(1,13)=Nnum(7);
            N(1,15)=Nnum(8);
            
            N(2,2)=Nnum(1);
            N(2,4)=Nnum(2);
            N(2,6)=Nnum(3);
            N(2,8)=Nnum(4);
            N(2,10)=Nnum(5);
            N(2,12)=Nnum(6);
            N(2,14)=Nnum(7);
            N(2,16)=Nnum(8);
            % Calcul du vecteur des forces volumiques
            % Initialisation
            f=zeros(2,1);
            % Contruction numérique du vecteur de sollicitation élémentaire
            % Terme des forces volumiques (poids, force magnétique, inertie,...)
            fe = fe  + W(i)*(N')*f*DJ;
        end
        % ------- Fin de la boucle sur les points d'intégration ------
        
        % Assemblage du ke en cours dans KG
        Kloc=Feeldof(nd,nnel,ddln); % Retourne les # de DDL concernés par l'élément
        KG=Assembl(KG,ke,Kloc);     % Assemble ke dans KG connaissant Kloc
        % Assemblage du vecteur sollicitations élémentaire fe
        FG= AssemblF(FG,fe,Kloc);   % Assemble fe dans FG connaissant Kloc
        Kloc2=Feeldof(nd,nnel,1);
        MG=Assembl(MG,me,Kloc2);     % Assemble me dans MG connaissant Kloc2 %matrice masse globale
    end
elseif dim=='Q9'
    for ie=1:nelin
        % Extraire les infos propres à cet élément (coord., connect.,propriétés. ...)
        % Boucle sur les noeuds de l'élément en question
        for i=1:nnel
            nd(i) = matConnec(ie,i);        % Connectivité des noeuds de l'élément
            coor(i,:) = matNoeud(nd(i),:);  % Coordonnées des noeuds de l'élément
        end
        coor;
        nd;
        % Initialisation de la matrice de rigidité élémentaire
        ke=zeros(Nui,Nui);
        me= zeros(nnel,nnel);
        % Initialisation du vecteur sollicitation élémentaire
        fe=zeros(Nui,1);
        % --------- Boucle sur les points d'intégration ----------
        for i=1:npg
            % Extraction des coordonnees des points d'intégration
            k = xloc(i,1);  % Ksi du point d'intégration sur le carré en cours
            n = xloc(i,2);  % Eta du point d'intégration sur le carré en cours
            % Évaluation numérique de la matrice Bksieta pour le point et l'élément
            Bke= BKsiEtaQ9(k,n);
            % Évaluation numérique des fonctions N pour le point et l'élément
            Nnum = InterpolQ9(k,n);
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
            % Calcul numérique de la matrice Bxy
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
            % Stockage de la matrice B au point d'intégration en cours
            % Note: cette matrice sera nécessaire pour le calcul des contraintes
            %    gradient(:,:,(ie-1)*(npg)+i)=B(:,:);
            % Contruction numérique de la matrice de rigidite elementaire
            % Note: Sommation des (poids) X (fonction évaluée au point de Gauss)
            ke = ke  + W(i)*(B')*D*B*DJ;
            me= me + W(i)*(Nnum')*Nnum*DJ;  % matrice masse qui va servir pour le lissage
            % Calcul de la matrice [N](2x18 dans le cas des Q9) à partir de Nnum
            % Initialisation
            N=zeros(2,nnel*2);
            % Assignation directe des 18 termes de N non nuls
            N(1,1)=Nnum(1);
            N(1,3)=Nnum(2);
            N(1,5)=Nnum(3);
            N(1,7)=Nnum(4);
            N(1,9)=Nnum(5);
            N(1,11)=Nnum(6);
            N(1,13)=Nnum(7);
            N(1,15)=Nnum(8);
            N(1,17)=Nnum(9);
            
            N(2,2)=Nnum(1);
            N(2,4)=Nnum(2);
            N(2,6)=Nnum(3);
            N(2,8)=Nnum(4);
            N(2,10)=Nnum(5);
            N(2,12)=Nnum(6);
            N(2,14)=Nnum(7);
            N(2,16)=Nnum(8);
            N(2,18)=Nnum(9);
            % Calcul du vecteur des forces volumiques
            % Initialisation
            f=zeros(2,1);
            % Contruction numérique du vecteur de sollicitation élémentaire
            % Terme des forces volumiques (poids, force magnétique, inertie,...)
            fe = fe  + W(i)*(N')*f*DJ;
        end
        % ------- Fin de la boucle sur les points d'intégration ------
        % Assemblage du ke en cours dans KG
        Kloc=Feeldof(nd,nnel,ddln); % Retourne les # de DDL concernés par l'élément
        KG=Assembl(KG,ke,Kloc);     % Assemble ke dans KG connaissant Kloc
        % Assemblage du vecteur sollicitations élémentaire fe
        FG= AssemblF(FG,fe,Kloc);   % Assemble fe dans FG connaissant Kloc
        Kloc2=Feeldof(nd,nnel,1); % Retourne les # de DDL concernés par l'élément
        MG=Assembl(MG,me,Kloc2);     % Assemble me dans MG connaissant Kloc2 %matrice masse globale
    end
elseif dim=='T6'
    for ie=1:nelin
    % Extraire les infos propres à cet élément (coord., connect.,propriétés. ...)
    % Boucle sur les noeuds de l'élément en question
    for i=1:nnel
        nd(i) = matConnec(ie,i);        % Connectivité des noeuds de l'élément
        coor(i,:) = matNoeud(nd(i),:);  % Coordonnées des noeuds de l'élément
    end
    coor;
    nd;
    % Initialisation de la matrice de rigidité élémentaire
    ke=zeros(Nui,Nui);
    me=zeros(nnel,nnel);
    %me= zeros(nnel,nnel);
    % Initialisation du vecteur sollicitation élémentaire
    fe=zeros(Nui,1);
    % --------- Boucle sur les points d'intégration ----------
    for i=1:npg
        % Extraction des coordonnees des points d'intégration
        k = xloc(1,i);  % Ksi du point d'intégration sur le carré en cours
        n = xloc(2,i);  % Eta du point d'intégration sur le carré en cours
        % Évaluation numérique de la matrice Bksieta pour le point et l'élément
        Bke= BKsiEtaT6(k,n);
        % Évaluation numérique des fonctions N pour le point et l'élément
        Nnum = InterpolT6(k,n);
        %xg=Nnum*coor(:,1);
        %yg=Nnum*coor(:,2);
        %plot(xg,yg,'.r')
        % Calcul numérique de la matrice Jacobienne (élément isoparamétrique)
        J= Bke*coor;
        % Calcul de son inverse
        j= J^-1 ;
        % Calcul de son déterminant
        %(utile car l'intégration se fait sur élément de référence plus loin)
        DJ= det(J);
        % Calcul numérique de la matrice Bxy
        Bxy = j*Bke;
        % Construction de la matrice B car élasticité 2D (avec les termes de Bxy)
        % RAPPEL: {epsilon}e=[B]*{U}e
        % [B](3*16) fait le lien entre les 16 déplacements et les 3 déformations
        % Initialisation
        B=zeros(3,nnel*2);
        % Assignation directe des 112 termes de B non nuls
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
        % Stockage de la matrice B au point d'intégration en cours
        % Note: cette matrice sera nécessaire pour le calcul des contraintes
        %    gradient(:,:,(ie-1)*(npg)+i)=B(:,:);
        % Contruction numérique de la matrice de rigidite elementaire
        % Note: Sommation des (poids) X (fonction évaluée au point de Gauss)
        ke = ke  + W(i)*(B')*D*B*DJ;
        me = me  + W(i)*(Nnum')*Nnum*DJ;
        % Calcul de la matrice [N](2x16 dans le cas des Q8) à partir de Nnum
        % Initialisation
        N=zeros(2,nnel*2);
        % Assignation directe des 16 termes de N non nuls
        N(1,1)=Nnum(1);
        N(1,3)=Nnum(2);
        N(1,5)=Nnum(3);
        N(1,7)=Nnum(4);
        N(1,9)=Nnum(5);
        N(1,11)=Nnum(6);
        
        N(2,2)=Nnum(1);
        N(2,4)=Nnum(2);
        N(2,6)=Nnum(3);
        N(2,8)=Nnum(4);   
        N(2,10)=Nnum(5);
        N(2,12)=Nnum(6);
        % Calcul du vecteur des forces volumiques
        % Initialisation      
        f=zeros(2,1);
        % Contruction numérique du vecteur de sollicitation élémentaire
        % Terme des forces volumiques (poids, force magnétique, inertie,...)
        fe = fe  + W(i)*(N')*f*DJ;
    end
    % ------- Fin de la boucle sur les points d'intégration ------
    
    % Assemblage du ke en cours dans KG
    Kloc=Feeldof(nd,nnel,ddln); % Retourne les # de DDL concernés par l'élément
    KG=Assembl(KG,ke,Kloc);     % Assemble ke dans KG connaissant Kloc
    % Assemblage du vecteur sollicitations élémentaire fe
    FG= AssemblF(FG,fe,Kloc);   % Assemble fe dans FG connaissant Kloc
    Kloc2=Feeldof(nd,nnel,1); % Retourne les # de DDL concernés par l'élément
    MG=Assembl(MG,me,Kloc2);     % Assemble me dans MG connaissant Kloc2 %matrice masse globale
    end
else
    disp('Element non reconnu')
end
end
