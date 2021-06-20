function Affiche_elements(cn,xy,L,h)
% Cette fonction affiche les �l�ments et leurs num�ros 
% Ouverture d'une nouvelle figure
figure
% Ne pas modifier les figures pr�c�dentes
hold on
% Calcul du nombre d'�l�ments
ne=length(cn(:,1));
% Boucle sur les �l�ments
for i=1:ne
    % D�limitation de l'aire de l'�l�ment
    xe=xy(cn(i,:),1);
    ye=xy(cn(i,:),2);
    % Remplir cette aire avec une couleur 
    fill(xe,ye,'g')
    % �crire au milieu de l'�l�ment son num�ro
    text(mean(xe),mean(ye),int2str(i))
end
% Ajustement du titre et noms des axes
title('�l�ments')
xlabel('x(m)')
ylabel('y(m)')
% Ajustement automatique des axes
axis([0 L -h h])
hold off
