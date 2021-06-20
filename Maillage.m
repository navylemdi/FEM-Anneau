function []= Maillage(xy,ne,nn,cn)

% xy=Gcoor;
% ne=nelt;
% nn=nnt;
% cn=Connec(1:nelt,:);

% ... affichage maillage(Gcoor,nelt,nnt,Connec(1:nelt,:));
if size(cn,2)~=9
    figure
    hold on
    axis('equal')
    title('maillage')
    plot(xy(:,1),xy(:,2),'.b','MarkerSize',20);              % ... grille des noeuds en .bleu
    
    % 1 - elements
    for i=1:ne
        hi=fill(xy(cn(i,:),1),xy(cn(i,:),2),'w');    % remplir zone white et relier
        set(hi,'LineWidth',1); % affichage contour elements par augmentation de la largeur
        set(hi,'FaceAlpha',0.1);
        xc=mean(xy(cn(i,:),1));              % localisation en x affichage no element (moy coor y)
        yc=mean(xy(cn(i,:),2));              % localisation en y affichage no element (moy coor y)
        si=text(xc,yc,int2str(i));           % localisation ... ? centrée surface
        set(si,'Color','r','FontSize',8);   % affichage no des elements
    end
    
    % 2 - noeuds
    for i=1:nn
        if xy(i,2)>=xy(end,2)/2
            si=text(xy(i,1),xy(i,2),[' ',int2str(i)],'VerticalAlignment','bottom');
        else
            si=text(xy(i,1),xy(i,2),[' ',int2str(i)],'VerticalAlignment','top');
        end
        set(si,'Color','b','FontSize',8);   % affichage no des noeuds
    end
    xlabel('x (m)')
    ylabel('y (m)')
    hold off
elseif size(cn,2)==9
    figure
    hold on
    axis('equal')
    title('maillage')
    plot(xy(:,1),xy(:,2),'.b','MarkerSize',20);              % ... grille des noeuds en .bleu
    
    % 1 - elements
    for i=1:ne
        hi=fill(xy(cn(i,1:8),1),xy(cn(i,1:8),2),'w');    % remplir zone white et relier
        set(hi,'LineWidth',2);               % affichage contour elements par augmentation de la largeur
        set(hi,'FaceAlpha',0.1);
        xc=mean(xy(cn(i,:),1));              % localisation en x affichage no element (moy coor y)
        yc=mean(xy(cn(i,:),2));              % localisation en y affichage no element (moy coor y)
        si=text(xc,yc,int2str(i));           % localisation ... ? centrée surface
        set(si,'Color','r','FontSize',8);   % affichage no des elements
    end
    
    % 2 - noeuds
    for i=1:nn
        if xy(i,2)>=xy(end,2)/2
            si=text(xy(i,1),xy(i,2),[' ',int2str(i)],'VerticalAlignment','bottom');
        else
            si=text(xy(i,1),xy(i,2),[' ',int2str(i)],'VerticalAlignment','top');
        end
        set(si,'Color','b','FontSize',8);   % affichage no des noeuds
    end
    xlabel('x (m)')
    ylabel('y (m)')
    hold off
else
    error('élément non pris en compte');
    
end
end
