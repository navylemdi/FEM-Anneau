function Affiche_resultats(cn,xy,uxy,amp,var,titre)
% Cette fonction permet d'afficher les r�sultats des calculs de
% d�placements ou de contraintes
% Ouverture d'une nouvelle figure
%if size(cn,2)==8 || size(cn,2)==6
% Ne pas modifier les figures pr�c�dentes
hold on 
% Poutre non-d�form�e
% Commande patch pilot�e selon: 'PropertyName',propertyvalue
%          Aires       Sommets       Donn�es               Couleur
temp=patch('faces',cn,'vertices',xy,'facevertexcdata',[],'facecolor','none');
% Lignes pointill�es pour la structure initiale
set(temp,'Linestyle',':');
% Poutre d�form�e (remplissage avec une �chelle de couleur)
%      Aires      Sommets                Donn�es               Couleur
patch('faces',cn,'vertices',xy+amp*uxy,'facevertexcdata',var,'facecolor','interp');
% Ajustement des couleurs de l'�chelle pour la grandeur physique
caxis([min(var),max(var)]);
% Affichage de l�gende de couleur � c�t� du graphique
colorbar;
% �criture du titre du graph
title(titre);
xlabel('x(m)');
ylabel('y(m)');
% Ajustement automatique des axes
axis('equal');
hold off;

% elseif size(cn,2)==9
%     lim=cn(1,9);
%     hold on 
% % Poutre non-d�form�e
% % Commande patch pilot�e selon: 'PropertyName',propertyvalue
% %          Aires       Sommets       Donn�es               Couleur
% temp=patch('faces',cn(:,1:8),'vertices',xy(:,:),'facevertexcdata',[],'facecolor','none');
% % Lignes pointill�es pour la structure initiale
% set(temp,'Linestyle',':');
% % Poutre d�form�e (remplissage avec une �chelle de couleur)
% %      Aires      Sommets                Donn�es               Couleur
% patch('faces',cn(:,1:8),'vertices',xy(:,:)+amp*uxy(:,:),'facevertexcdata',var(:),'facecolor','interp');
% % Ajustement des couleurs de l'�chelle pour la grandeur physique
% caxis([min(var),max(var)]);
% % Affichage de l�gende de couleur � c�t� du graphique
% colorbar;
% % �criture du titre du graph
% title(titre);
% xlabel('x(m)');
% ylabel('y(m)');
% % Ajustement automatique des axes
% axis('equal');
% hold off;
% else
%     error('�l�ment non pris en compte');
%     
% end
% end
