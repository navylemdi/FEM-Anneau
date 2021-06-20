function Affiche_elements(cn,xy,L,h)
% Cette fonction affiche les éléments et leurs numéros 
% Ouverture d'une nouvelle figure
figure
% Ne pas modifier les figures précédentes
hold on
% Calcul du nombre d'éléments
ne=length(cn(:,1));
% Boucle sur les éléments
for i=1:ne
    % Délimitation de l'aire de l'élément
    xe=xy(cn(i,:),1);
    ye=xy(cn(i,:),2);
    % Remplir cette aire avec une couleur 
    fill(xe,ye,'g')
    % Écrire au milieu de l'élément son numéro
    text(mean(xe),mean(ye),int2str(i))
end
% Ajustement du titre et noms des axes
title('Éléments')
xlabel('x(m)')
ylabel('y(m)')
% Ajustement automatique des axes
axis([0 L -h h])
hold off
