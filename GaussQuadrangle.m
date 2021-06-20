function [xloc,W] =GaussQuadrangle(npg)
% %Points et poids de Gauss pour un carré
w1(1)= 8/9;
w1(2)= 5/9;
w1(3)= 5/9;
x1(1)=0;
x1(2)=-sqrt(3/5);
x1(3)= sqrt(3/5);

if npg==1
    xloc(1,1)=0;
    xloc(1,2)=0;
    W(1)= 4;
end

if npg==4
    % Méthode produit à 2 X 2 points
    % 4 points intègrent exactement jusqu'à l'ordre 3
    W(1)= 1;
    W(2)= 1;
    W(3)= 1;
    W(4)= 1;
    xloc(1,1)=  1/sqrt(3);
    xloc(1,2)=  1/sqrt(3);
    xloc(2,1)=  -1/sqrt(3);
    xloc(2,2)=  1/sqrt(3);
    xloc(3,1)=  -1/sqrt(3);
    xloc(3,2)=  -1/sqrt(3);
    xloc(4,1)=  1/sqrt(3);
    xloc(4,2)=  -1/sqrt(3);
end
if npg==9
    % Méthode à 9 points intégrant jusqu'à l'ordre m=5 dans chaque direction
    k=0;
    for i=1:3
        for j=1:3
            k=k+1;
            W(k)=w1(i)*w1(j);
            xloc(k,1)=x1(i);
            xloc(k,2)=x1(j);
        end
    end
end
end