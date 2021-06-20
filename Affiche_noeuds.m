function Affiche_noeuds(xy,L,h)
% Cette fonction affiche les noeuds et leurs numéros
% Ouverture d'une nouvelle fenêtre de visualisation
figure
% Dessine un point vis-à-vis chaque noeud
plot(xy(:,1),xy(:,2),'.')
% Calcul du nombre de noeud
n=length(xy(:,1));
% Écriture du numéro de noeud vis-à-vis chaque noeud
for i=1:n % Pour tous les noeuds
    %       x       y           #de noeud
    text( xy(i,1), xy(i,2) , [' ',int2str(i)] )
end
% Ajuster le titre
title('Noeuds')
xlabel('x(m)')
ylabel('y(m)')
% Ajustement automatique des axes
axis([0 L -h h])


