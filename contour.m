function FG=contour(npg, xlocC, WC, nelin, nelc, nnelc, ddlnc, Nuic, matConnec, matNoeud,q,FG)
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
    coor;
    nd;
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
end
