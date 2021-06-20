function val_Q8 = interpQ4_to_Q8(val_Q4,params, dims)
    type_maille = params{1};
    Rmin = params{2};
    Rmax = params{3};
    T = params{4};
    
    nRadius = dims{1};
    nTheta = dims{2};
    
    assert(size(val_Q4,1)==((nRadius+1)*(nTheta+1)));
    % ===========================================
    % Interpolation de Q4 à Q8 des contraintes
    % ===========================================
    t_1 = linspace(0,T,nTheta+1);
    r_1 = linspace(Rmin,Rmax,nRadius+1);
    [R_mesh,T_mesh] = meshgrid(r_1,t_1);
    val_Q4_mesh = reshape(val_Q4,[nRadius+1,nTheta+1]);
    % Demi points horizontaux
    nb_pt = (nRadius+1)*(nTheta+1);
    nb_pt_2 = nb_pt + (nRadius+1)*(nTheta);
    tmin_2 = T*(1/(2*nTheta));
    tmax_2 = T*(1-(1/(2*nTheta)));
    rmin_2 = Rmin;
    rmax_2 = Rmax;
    t_2 = linspace(tmin_2,tmax_2,nTheta);
    r_2 = linspace(rmin_2,rmax_2,nRadius+1);
    [R_mesh_2,T_mesh_2] = meshgrid(r_2,t_2);
    % interpolation sur le bloc 2
    val_Q8_mesh_2 = interp2(R_mesh,T_mesh,val_Q4_mesh,R_mesh_2,T_mesh_2,'cubic');
    val_Q8_2 = reshape(val_Q8_mesh_2,[(nRadius+1)*(nTheta),1]);
    % Demi points verticaux
    nb_pt_3 = nb_pt_2 + (nRadius)*(nTheta+1);
    tmin_3 = 0;
    tmax_3 = T;
    rmin_3 = Rmin+(Rmax-Rmin)*(1/(2*nRadius));
    rmax_3 = Rmax-(Rmax-Rmin)*(1/(2*nRadius));
    t_3 = linspace(tmin_3,tmax_3,nTheta+1);
    r_3 = linspace(rmin_3,rmax_3,nRadius);
    [R_mesh_3,T_mesh_3] = meshgrid(r_3,t_3);
    % interpolation sur le bloc 3
    val_Q8_mesh_3 = interp2(R_mesh,T_mesh,val_Q4_mesh,R_mesh_3,T_mesh_3,'cubic');
    val_Q8_3 = reshape(val_Q8_mesh_3,[(nRadius)*(nTheta+1),1]);
    % Assemblage du vecteur de sortie
    val_Q8 = zeros(nb_pt_3,1);
    val_Q8(1:nb_pt,1)=val_Q4;
    val_Q8(nb_pt+1:nb_pt_2,1)=val_Q8_2;
    val_Q8(nb_pt_2+1:nb_pt_3,1)=val_Q8_3;
end