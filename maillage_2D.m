function [Coord, Connect, ContourB, ContourH, ContourG, ContourD] = maillage_2D(params, dims)
    
    % vérification des arguments
    assert(length(params)==4)
    assert(length(dims)==2)
    
    type_maille = params{1};
    Rmin = params{2};
    Rmax = params{3};
    T = params{4};
    
    nRadius = dims{1};
    nTheta = dims{2};
    
    % Coordonnées
    Coord = construct_coord(nRadius,nTheta,T,Rmin,Rmax,type_maille);  
     
    % Matrice de connexion
    Connect = construct_connect(nRadius,nTheta,type_maille);  
    % Contour
    [cntr_b,cntr_g,cntr_h,cntr_d] = construct_contour(nRadius,nTheta);
    ContourB=zeros(nTheta,size(Connect,2));
    ContourH=zeros(nTheta,size(Connect,2));
    ContourG=zeros(nRadius,size(Connect,2));
    ContourD=zeros(nRadius,size(Connect,2));
    l=1;
    for k=1:nTheta
        ContourB(k,1:3)=cntr_b(l:l+2);%permet de remette bien les valeurs 
        ContourH(k,1:3)=cntr_h(l:l+2);
        l=l+2;
    end   
    l=1;
    for k=1:nRadius
        ContourG(k,1:3)=cntr_g(l:l+2);%permet de remette bien les valeurs 
        ContourD(k,1:3)=cntr_d(l:l+2);
        l=l+2;
    end 
    %ContourG = [cntr_g',zeros(length(cntr_g),1)];
    %ContourD = [cntr_d',zeros(length(cntr_d),1)];
    [Coord(:,1), Coord(:,2)] = pol2cart(Coord(:,2),Coord(:,1));
    Coord=[Coord(:,1), Coord(:,2)];
    
    function Coord = construct_coord(nRadius,nTheta,T,Rmin,Rmax,type_maille)
        % Coins des éléments
        t = linspace(0,T,nTheta+1);
        r = linspace(Rmin,Rmax,nRadius+1);
        [R_rshp,T_rshp] = meshgrid(r,t);
        nb_pt = (nTheta+1)*(nRadius+1);
        Coord(1:nb_pt,1) = reshape(R_rshp,[nb_pt,1]);
        Coord(1:nb_pt,2) = reshape(T_rshp,[nb_pt,1]);
        % Demi points horizontaux
        nb_pt_2 = nb_pt + (nRadius+1)*(nTheta);
        tmin_2 = T*(1/(2*nTheta));
        tmax_2 = T*(1-(1/(2*nTheta)));
        t_2 = linspace(tmin_2,tmax_2,nTheta);
        rmin_2 = Rmin;
        rmax_2 = Rmax;
        r_2 = linspace(rmin_2,rmax_2,nRadius+1);
        [R_rshp_2,T_rshp_2] = meshgrid(r_2,t_2);
        Coord(nb_pt+1:nb_pt_2,1) = reshape(R_rshp_2,[(nRadius+1)*(nTheta),1]);
        Coord(nb_pt+1:nb_pt_2,2) = reshape(T_rshp_2,[(nRadius+1)*(nTheta),1]);
        % Demi points verticaux
        nb_pt_3 = nb_pt_2 + (nRadius)*(nTheta+1);
        tmin_3 = 0;
        tmax_3 = T;
        t_3 = linspace(tmin_3,tmax_3,nTheta+1);
        rmin_3 = Rmin+(Rmax-Rmin)*(1/(2*nRadius));
        rmax_3 = Rmax-(Rmax-Rmin)*(1/(2*nRadius));
        r_3 = linspace(rmin_3,rmax_3,nRadius);
        [R_rshp_3,T_rshp_3] = meshgrid(r_3,t_3);
        Coord(nb_pt_2+1:nb_pt_3,1) = reshape(R_rshp_3,[(nRadius)*(nTheta+1),1]);
        Coord(nb_pt_2+1:nb_pt_3,2) = reshape(T_rshp_3,[(nRadius)*(nTheta+1),1]); 
        % Points centraux pour les éléments Q9 et T6
        if strcmp(type_maille,'Q9') || strcmp(type_maille,'T6')
            nb_pt_add = nb_pt_3 + (nRadius)*(nTheta);
            tmin_add = T*(1/(2*nTheta));
            tmax_add = T*(1-(1/(2*nTheta)));
            t_add = linspace(tmin_add,tmax_add,nTheta);
            Rmin_add = Rmin+(Rmax-Rmin)*(1/(2*nRadius));
            Rmax_add = Rmax-(Rmax-Rmin)*(1/(2*nRadius));
            r_add = linspace(Rmin_add,Rmax_add,nRadius);
            [R_rshp_add,T_rshp_add] = meshgrid(r_add,t_add);
            Coord(nb_pt_3+1:nb_pt_add,1) = reshape(R_rshp_add,[(nRadius)*(nTheta),1]);
            Coord(nb_pt_3+1:nb_pt_add,2) = reshape(T_rshp_add,[(nRadius)*(nTheta),1]);
        end
        % Affichage des noeuds pour le debug
        affiche = 0;
        if affiche
            plot(Coord(1:nb_pt,1),Coord(1:nb_pt,2),'ro')
            hold on
            plot(Coord(nb_pt+1:nb_pt_2,1),Coord(nb_pt+1:nb_pt_2,2),'k|')
            plot(Coord(nb_pt_2+1:nb_pt_3,1),Coord(nb_pt_2+1:nb_pt_3,2),'k_')
            if strcmp(type_maille,'Q9') || strcmp(type_maille,'T6')
                plot(Coord(nb_pt_3+1:nb_pt_add,1),Coord(nb_pt_3+1:nb_pt_add,2),'kd')
            end
        end
    end
    
    function [cntr_b,cntr_g,cntr_h,cntr_d] = construct_contour(nRadius,nTheta)
        nb_case = nRadius*nTheta;
        nb_coin = (nTheta+1)*(nRadius+1);
        nb_demih = (nRadius+1)*(nTheta)+nb_coin;

        cntr_b = zeros(1,2*nTheta+1);
        cntr_b = zeros(1,2*nTheta+1);
        cntr_g = zeros(1,2*nRadius+1);
        cntr_d = zeros(1,2*nRadius+1);

        for i =1:nTheta
            % Contour du rayon intérieur
            cntr_b(2*i-1) = i;
            cntr_b(2*i) = nb_coin+i;
            % Contour du rayon extérieur
            cntr_h(2*i-1) = nb_coin-(nTheta+1)+i;
            cntr_h(2*i) = nb_demih-nTheta+i;
        end
        cntr_b(2*nTheta+1) = nTheta+1;
        cntr_h(2*nTheta+1) = nb_coin;
        for j = 1:nRadius
            % Contour gauche
            cntr_g(2*j-1) = 1+(j-1)*(nTheta+1);
            cntr_g(2*j) = nb_demih+(j-1)*(nTheta+1)+1;
            % Contour droite
            cntr_d(2*j-1) = j*(nTheta+1);
            cntr_d(2*j) = nb_demih+j*(nTheta+1);
        end
        cntr_g(2*nRadius+1) = nb_coin-nTheta;
        cntr_d(2*nRadius+1) = nb_coin;
    end
    
    function Connect = construct_connect(nRadius,nTheta,type_maille)
        nb_case = nRadius*nTheta;
        nb_coin = (nTheta+1)*(nRadius+1);
        nb_demih = (nRadius+1)*(nTheta)+nb_coin;
        nb_demiv = (nRadius)*(nTheta+1)+nb_demih;

        if strcmp(type_maille,'T6')
            Connect = zeros(2*nb_case,6);
            for j =1:nRadius
                for i = 1:nTheta
                    ind_case = i + (j-1)*nTheta ;
                    ind_pt = i+(j-1)*(nTheta+1);
                    if i > nTheta/2
                        %élément T6 impair
                        Connect(2*ind_case-1,1) = ind_pt;
                        Connect(2*ind_case-1,2) = ind_pt+1;
                        Connect(2*ind_case-1,3) = ind_pt+(nTheta+1);
                        Connect(2*ind_case-1,4) = nb_coin+ind_case;
                        Connect(2*ind_case-1,5) = nb_demiv+ind_case;
                        Connect(2*ind_case-1,6) = nb_demih+ind_pt;
                        %élément T6 pair
                        Connect(2*ind_case,1) = ind_pt+1;
                        Connect(2*ind_case,2) = ind_pt+(nTheta+1)+1;
                        Connect(2*ind_case,3) = ind_pt+(nTheta+1);
                        Connect(2*ind_case,4) = nb_demih+ind_pt+1;
                        Connect(2*ind_case,5) = nb_coin+ind_case+(nTheta);
                        Connect(2*ind_case,6) = nb_demiv+ind_case;
                    else
                        %élément T6 impair
                        Connect(2*ind_case-1,1) = ind_pt;
                        Connect(2*ind_case-1,2) = ind_pt+1;
                        Connect(2*ind_case-1,3) = ind_pt+(nTheta+1)+1;
                        Connect(2*ind_case-1,4) = nb_coin+ind_case;
                        Connect(2*ind_case-1,5) = nb_demih+ind_pt+1;
                        Connect(2*ind_case-1,6) = nb_demiv+ind_case;
                        
                        %élément T6 pair
                        
                        Connect(2*ind_case,1) = ind_pt+(nTheta+1)+1;
                        Connect(2*ind_case,2) = ind_pt+(nTheta+1);
                        Connect(2*ind_case,3) = ind_pt;
                        
                        Connect(2*ind_case,4) = nb_coin+ind_case+(nTheta);
                        Connect(2*ind_case,5) = nb_demih+ind_pt;
                        Connect(2*ind_case,6) = nb_demiv+ind_case;
                    end
                end
            end
        elseif strcmp(type_maille,'Q8')
            Connect = zeros(nb_case,8);
            for j =1:nRadius
                for i = 1:nTheta
                    ind_case = i + (j-1)*nTheta ;
                    ind_pt = i+(j-1)*(nTheta+1);
                    % Élément Q8
                    Connect(ind_case,1) = ind_pt;
                    Connect(ind_case,2) = ind_pt+1;
                    Connect(ind_case,3) = ind_pt+(nTheta+1)+1;
                    Connect(ind_case,4) = ind_pt+(nTheta+1);
                    Connect(ind_case,5) = nb_coin+ind_case;
                    Connect(ind_case,6) = nb_demih+ind_pt+1;
                    Connect(ind_case,7) = nb_coin+ind_case+(nTheta);
                    Connect(ind_case,8) = nb_demih+ind_pt;

                end
            end
        elseif strcmp(type_maille,'Q9')
            Connect = zeros(nb_case,9);
            for j =1:nRadius
                for i = 1:nTheta
                    ind_case = i + (j-1)*nTheta ;
                    ind_pt = i+(j-1)*(nTheta+1);
                    % Élément Q9
                    Connect(ind_case,1) = ind_pt;
                    Connect(ind_case,2) = ind_pt+1;
                    Connect(ind_case,3) = ind_pt+(nTheta+1)+1;
                    Connect(ind_case,4) = ind_pt+(nTheta+1);
                    Connect(ind_case,5) = nb_coin+ind_case;
                    Connect(ind_case,6) = nb_demih+ind_pt+1;
                    Connect(ind_case,7) = nb_coin+ind_case+(nTheta);
                    Connect(ind_case,8) = nb_demih+ind_pt;
                    Connect(ind_case,9) = nb_demiv+ind_case;
                end
            end
        else
            error('éléments non pris en compte')
        end
    end
end