function [uv,uv_th, tab_sigma, tab_sigma_th, Coord, Connect] = modele_Q9(nRadius,nTheta)
    % Structure en équilibre statique
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ********** PARAMÈTRES DU MODÈLE **********
    % Propriétés du matériau
    E = 7e10;  % Module de Young en Pa
    nu = 0.3 ; % Coefficient de Poisson
    e= 1e-3;  % épaisseur en m selon z
    a = 10e-3; % rayon interne m 
    b = 20e-3; % rayon externe m 

    % ********* CHARGEMENT DE LA POUTRE *********
    % Charge répartie
    q= 1e9;%pression radiale N/m2

    % ********** Paramètre maillage *********
    % Lecture du nombre de ddl par noeud
    ddln = 2;
    % Lecture du nombre de ddl par noeud sur les contours
    ddlnc = 2;
    % Lecture du nombre de noeud par élément interne
    nnel = 9;%6%8%9 en fonction
    % Nombre de noeud de l'élément de contour
    nnelc=3;
    % Calcul du nombre de DDL par élément interne (18 dans le cas du Q9)
    Nui=ddln*nnel;
    % Calcul du nombre de DDL par élément de contour (6 pour barre quadratique)
    Nuic=ddlnc*nnelc;
    % filename='xy1x3Q8.dat';
    % dim=2;
    % typeAnal='structure';
    %Dimensionnement du maillage
    type_maille = 'Q9';
    %Paramètre du problème
    T = pi/2;
    Rmin = a;
    Rmax = b;

    params = {type_maille,Rmin,Rmax,T};
    dims = {nRadius,nTheta};
    %[matNoeud, matConnec, DirichletX ,DirichletY, RotZ] = ansys2matlab(filename,dim,typeAnal);
    [matNoeud, matConnectelemint, ContourB, ContourH, DirichletY,DirichletX] = maillage_2D(params, dims);

    matConnec=zeros(size(matConnectelemint,1)+size(ContourB,1),size(matConnectelemint,2));
    matConnec(1:size(matConnectelemint,1),:)=matConnectelemint;
    matConnec(size(matConnectelemint,1)+1:end,:)=ContourB;

    nelin=0;
    for i=1:size(matConnec,1)
        if (matConnec(i,4)~=0)
            nelin=nelin+1;
        end
    end
    nelc=size(matConnec,1)-nelin;%Nombre total d'éléments du contour
    nnt=length(matNoeud(:,1));%Nombre de noeuds total
    ndlt= nnt*ddln;%Nombre total de d.d.l
    %On mets les connectivités dans le bon ordre
    NewmatConnec=zeros(nelin,nnel-1);
    NewmatConnec(1:nelin,1:2:8)=matConnec(1:nelin,1:4);
    NewmatConnec(1:nelin,2:2:8)=matConnec(1:nelin,4+1:8);
    matConnec(1:nelin,1:8)=NewmatConnec;
    %Maillage(matNoeud,nelin,nnt,matConnec);
    for i=1:nelin
        for j=1:4
            c=matConnec(i,j);
            d=matConnec(i,8-j+1);
            matConnec(i,8-j+1)=c;
            matConnec(i,j)=d;
        end
    end
    for i=nelin+1:nelin+nelc
        for j=1:3
            c=matConnec(i,1);
            d=matConnec(i,3);
            matConnec(i,3)=c;
            matConnec(i,1)=d;
        end
    end
    matConnec(1:nelin,1:8)=circshift(matConnec(1:nelin,1:8),1,2);

    %====================================================================
    % 4. Initialisation de la matrice de rigidité et du vecteur f globaux
    %====================================================================
    KG= sparse(ndlt,ndlt);
    MG= sparse(nnt,nnt);  % matrice masse globale

    FG= zeros(ndlt,1);
    %=============================================================
    % 5. Matrice des propriétés physiques CONTRAINTE PLANE Sz = 0
    %=============================================================
    %if(Hypothesis  ==1)    %plane stress
    facteur=E/(1-nu^2);
    D=facteur*[ 1 , nu , 0 ;
                nu , 1 , 0 ;
                0 , 0 , (1-nu)/2];
    %=============================================================
    % 6.        Calculs préliminaires à la boucle
    %=============================================================
    % Points d'intégration de Gauss
    % Utilisation de l'intégration numérique 3 X 3 symétrique sur les
    % quadrilataires
    npg = 9;
    [xloc,W]= GaussQuadrangle(npg);%Ou Quadrangle
    %=============================================================
    %************ 7. Boucle sur les éléments internes ************
    %=============================================================

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
    % *****Fin de la boucle sur les éléments intérieurs*****
    %=============================================================
    %*********  8. Boucle sur les éléments de contour ************
    %=============================================================
    % Ajustement des paramètres d'intégration
    npg=4;  % Nombre de points de Gauss sur les éléments barres
    [xlocC,WC]= Gauss(npg);%Ou Quadrangle
    % Boucle sur les barres de contour
    clear nd coor
    for ie=1+nelin:nelc+nelin
        % (Car on a mis leurs connectivités à la suite des Q8  dans Connec)
        % Initialiser à zéro le vecteur élémentaire fe
        fe=zeros(Nuic,1);
        % Extraire les infos propres à cet élément (coord., connect.,propriétés. ...)
        % Boucle sur les noeuds de l'élément en question
        for i=1:nnelc
            nd(i) = matConnec(ie,i);% Connectivités des noeuds de l'élément
            coor(i,:) = matNoeud(nd(i),:);% Coordonnées des noeuds de l'élément
        end
        % Longueur de l'élément par le théorème de Pythagore
        Le=  (coor(1,1)-coor(2,1))^2+ (coor(1,2)-coor(2,2))^2;
        Le = 2*sqrt(Le);
        % Calcul du déterminant de la matrice J
        DJ=Le/2;
        %------- Boucle sur les points d'intégration -----------
        for i=1:npg
            % Extraction des coordonnées d'intégration
            k = xlocC(i,1); % Ksi du point d'intégration sur la barre en cours
            % Évaluation numérique des fonctions N pour le point et l'élément
            Nnum = Interpolc(k);
            Bke=BKsiEtac(k);
            %on recupere la position réelle des points d'intégration pour
            %calculer les forces composée dû à la pression
            xg=Nnum*coor(1:3,1);%car on change que le 3 premieres lignes de coor
            yg=Nnum*coor(1:3,2);
            theta=atan(yg/xg);
            %theta*(180/pi);
            % Calcul de la matrice [N](2x6 dans le cas des barres) avec Nnum
            % Initialisation
            N=zeros(2,nnelc*2);
            % Assignation directe des 6 termes de N non nuls
            N(1,1)=Nnum(1);
            N(1,3)=Nnum(2);
            N(1,5)=Nnum(3);

            N(2,2)=Nnum(1);
            N(2,4)=Nnum(2);
            N(2,6)=Nnum(3);
            % Calcul du vecteur des forces surfaciques
            % Initialisation
            fs=zeros(2,1);
            % Calculs des termes du vecteur
            fs(1)=q*cos(theta);
            fs(2)=q*sin(theta);
            % Contruction numérique du vecteur de charge élémentaire
            % Il s'agit ici des forces surfaciques (charge répartie, etc.)
            fe = fe  + WC(i)*(N')*fs*DJ;
        end
        % ------ Fin de la boucle sur les points d'intégration ---
        % Assemblage dans le système global aux endroits appropriés
        Kloc=Feeldof(nd,nnelc,ddlnc);
        FG= AssemblF(FG,fe,Kloc);
    end
    %=============================================================
    % 9. Traitement des C.L. de déplacement imposé (Dirichlet)
    %=============================================================
    % Imposition des déplacements u (selon x)
    nclx = size(DirichletX,1);
    for i=1:nclx          % Boucle sur les noeuds avec u imposé
        for j=1:3
            nn= DirichletX(i,j);     % Numéro du noeud en question
            % Calcul de la position de ce DDL
            pos=(nn-1)*ddln+1;  % u et v des noeuds avant + 1
            for l=1:ndlt        % Boucle sur tous les DDL du système
                if l==pos
                    KG(pos,l)= 1 ;     % Met un 1 si c'est le bon DDL
                    FG(pos)= 0;  % Impose u dans le terme de droite
                else
                    KG(pos,l)=0;    % Met un 0 si ce n'est pas le bon DDL
                end
            end
        end
    end

    % Imposition des déplacements v (selon y)
    ncly = size(DirichletY,1);
    for i=1:ncly            % Boucle sur les noeuds avec v imposé
        for j=1:3
            nn= DirichletY(i,j);       % Numéro du noeud en question
            % Calcul de la position de ce DDL
            pos=(nn-1)*ddln+2;  % u et v des noeuds avant + 2
            for l=1:ndlt        % Boucle sur tous les DDL du système
                if l==pos
                    KG(pos,l)= 1 ;     % Met un 1 si c'est le bon DDL
                    FG(pos)= 0;  % Impose u dans le terme de droite
                else
                    KG(pos,l)=0;    % Met un 0 si ce n'est pas le bon DDL
                end
            end
        end
    end

    sol= KG\FG;
    cn=matConnec(1:nelin,1:8);  % Connectivités des quadrilatères seulement
    u=sol(1:2:2*nnt); % Déplacements en x en m
    v=sol(2:2:2*nnt);     % Déplacements en y en m

    %===============================================================
    % 11. CALCUL DES CONTRAINTES Aux point de Gauss puis on fait un lissage pour les obtenir aux noeuds
    %===============================================================
    GX= zeros(nnt,1);  % vecteur auxiliare
    GY= zeros(nnt,1);
    GXY= zeros(nnt,1);
    % Ajustement du npg pour les éléments internes
    npg=9;
    [xloc,W]= GaussQuadrangle(npg);
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

    % Obtenir les contraintes lissees aux noeuds
    sigmax= MG\GX;
    sigmay= MG\GY;
    sigmaxy= MG\GXY;
    contraintes_n(1,:)=sigmax';
    contraintes_n(2,:)=sigmay';
    contraintes_n(3,:)=sigmaxy';

    sigmar=zeros(length(sigmax),1);
    sigmat=zeros(length(sigmay),1);
    sigmart=zeros(length(sigmaxy),1);
    sigmarth=zeros(length(sigmax),1);
    sigmatth=zeros(length(sigmay),1);
    sigmartth=zeros(length(sigmaxy),1);
    sigmaxth=zeros(length(sigmax),1);
    sigmayth=zeros(length(sigmay),1);
    sigmaxyth=zeros(length(sigmaxy),1);
    for i=1:size(matNoeud,1)
        x=matNoeud(i,1);
        y=matNoeud(i,2);
        theta=atan(y/x);
        theta*180/pi;
        R=sqrt(x^2+y^2);
        sigmar(i)=sigmax(i)*cos(theta)^2+sigmay(i)*sin(theta)^2+sigmaxy(i)*sin(2*theta);
        sigmat(i)=sigmax(i)*sin(theta)^2+sigmay(i)*cos(theta)^2-sigmaxy(i)*sin(2*theta);
        sigmart(i)=cos(theta)*sin(theta)*(sigmay(i)-sigmax(i))+sigmaxy(i)*cos(2*theta);
        sigmarth(i)=a^2*(1-b^2/R^2)*q/(b^2-a^2);
        sigmatth(i)=a^2*(1+b^2/R^2)*q/(b^2-a^2);
        sigmartth(i)=0;
        sigmaxth(i)=sigmarth(i)*cos(theta)^2+sigmatth(i)*sin(theta)^2+sigmartth(i)*sin(2*theta);
        sigmayth(i)=sigmarth(i)*sin(theta)^2+sigmatth(i)*cos(theta)^2-sigmartth(i)*sin(2*theta);
        sigmaxyth(i)=-(cos(theta)*sin(theta)*(sigmatth(i)-sigmarth(i))+sigmartth(i)*cos(2*theta));
        
        Y = b/a;
        g_th(i,1) = ((1-nu)/E)*q*R/(Y^2-1) + ((1+nu)/E)*q*b^2/((Y^2-1)*R);
        u_th(i,1) = g_th(i,1)*cos(theta);
        v_th(i,1) = g_th(i,1)*sin(theta);
    end
    

    
    
    
    uv = [u,v];
    uv_th = [u_th,v_th];
    tab_sigma = [sigmax,sigmay,sigmaxy,sigmar,sigmat,sigmart];
    tab_sigma_th = [sigmaxth,sigmayth,sigmaxyth,sigmarth,sigmatth,sigmartth];
    Coord = matNoeud;
    Connect = matConnec;
end