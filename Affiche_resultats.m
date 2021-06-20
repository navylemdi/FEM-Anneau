function Affiche_resultats(cn,xy,uxy,amp,var,titre)
% Cette fonction permet d'afficher les résultats des calculs de
% déplacements ou de contraintes
% Ouverture d'une nouvelle figure
%if size(cn,2)==8 || size(cn,2)==6
% Ne pas modifier les figures précédentes
hold on 
% Poutre non-déformée
% Commande patch pilotée selon: 'PropertyName',propertyvalue
%          Aires       Sommets       Données               Couleur
temp=patch('faces',cn,'vertices',xy,'facevertexcdata',[],'facecolor','none');
% Lignes pointillées pour la structure initiale
set(temp,'Linestyle',':');
% Poutre déformée (remplissage avec une échelle de couleur)
%      Aires      Sommets                Données               Couleur
patch('faces',cn,'vertices',xy+amp*uxy,'facevertexcdata',var,'facecolor','interp');
% Ajustement des couleurs de l'échelle pour la grandeur physique
caxis([min(var),max(var)]);
% Affichage de légende de couleur à côté du graphique
colorbar;
% Écriture du titre du graph
title(titre);
xlabel('x(m)');
ylabel('y(m)');
% Ajustement automatique des axes
axis('equal');
hold off;

% elseif size(cn,2)==9
%     lim=cn(1,9);
%     hold on 
% % Poutre non-déformée
% % Commande patch pilotée selon: 'PropertyName',propertyvalue
% %          Aires       Sommets       Données               Couleur
% temp=patch('faces',cn(:,1:8),'vertices',xy(:,:),'facevertexcdata',[],'facecolor','none');
% % Lignes pointillées pour la structure initiale
% set(temp,'Linestyle',':');
% % Poutre déformée (remplissage avec une échelle de couleur)
% %      Aires      Sommets                Données               Couleur
% patch('faces',cn(:,1:8),'vertices',xy(:,:)+amp*uxy(:,:),'facevertexcdata',var(:),'facecolor','interp');
% % Ajustement des couleurs de l'échelle pour la grandeur physique
% caxis([min(var),max(var)]);
% % Affichage de légende de couleur à côté du graphique
% colorbar;
% % Écriture du titre du graph
% title(titre);
% xlabel('x(m)');
% ylabel('y(m)');
% % Ajustement automatique des axes
% axis('equal');
% hold off;
% else
%     error('élément non pris en compte');
%     
% end
% end
