%%Interface ANSYS Mechanical vers format std Matlab
function [matNoeud, matConnec, DirichletX,DirichletY,DirichletZ] = ansys2matlab(fileName,dim,typeAnal)
%Coord: les coordonnes pour chacun des noeuds, nbNoeud x dim
%Connec: les connectivite des elements, ensuite les elements de contour,
%nbElement+nbElementContour x nbNoeudMax
%Dirichlet[X-Y-Z: les point ou l'on impose des conditions de dirichlet, nbPoints x [valeurX, valeurY]
%filename : le nom du fichier ds.dat dans le dossier MECH d'ANSYS
%dim le nomdre de dimension du problème 2 ou 3
%%
% 1.Presentement seulement les conditions limites de deplacement imposé et de
% force sur une surface est considéré
% 2. Le code doit être validé pour un maillage en 3dimension
% 3. Pour avoir les éléments et les éléments de contour de type linéaire,
% il faut enlever dans votre projet ANSYS Mechanical->analyse->... enlever
% les K:\DEVOIR4\ANSYS\DEVOIR4_files\dp0\SYS\MECH
% point milieux
%%
if ( size(findstr(typeAnal, 'therm')) ~= [0,0] )
    fprintf('Le fichier lu est pour une analyse thermique.\n')
else
    fprintf('Le fichier lu est pour une analyse de structure.\n')
end
%pause
fid = fopen(fileName,'r');
tline = fgetl(fid);
%On lit le fichier jusqu'a la fin
while ischar(tline)
    %%Pour lire les noeuds
    if (size(findstr(tline, 'nblock')) ~= [0,0])
        %Sautons l'entete
        tline = fgetl(fid);
        %lisons la premiere ligne
        tline = fgetl(fid);
        while ( size(tline,2) ~= 2 || size(findstr(tline, '-1'),1) == 0)
            tempNoeud = sscanf(tline,'%f');
            matNoeud(tempNoeud(1),:) = tempNoeud(2:1+dim);
            tline = fgetl(fid);
        end
        tline = fgetl(fid);
    end
    %%Pour lire les connections des éléments
    if (size(findstr(tline, 'eblock')) ~= [0,0])
        %Sautons l'entete
        tline = fgetl(fid);
        %lisons la premiere ligne
        tline = fgetl(fid);
        
        while ( size(tline,2) ~= 2 || size(findstr(tline, '-1'),1) == 0)
            
            tempConnec = sscanf(tline,'%f');
            
            if(dim == 2)%elements en 2 dimensions
                %pour les quad lin en
                if(tempConnec(9) == 4 && tempConnec(14) ~= tempConnec(15))
                    matConnec(tempConnec(11),:) = tempConnec(12:15);
                end
                %pour les quad quadratique incomplet
                if(tempConnec(9) == 8)
                    matConnec(tempConnec(11),:) = tempConnec(12:19);
                end
                %pour les triangle lin
                if(tempConnec(9) == 4 && tempConnec(14) == tempConnec(15) )
                    matConnec(tempConnec(11),:) = tempConnec(12:14);
                end
                %pour les triangle quadratique
                if(tempConnec(9) == 6)
                    matConnec(tempConnec(11),:) = tempConnec(12:17);
                end
            elseif (dim ==3)% elements en 3 dimensions
                %pour les hexa lin
                if(tempConnec(9) == 8)
                    matConnec(tempConnec(11),:) = tempConnec(12:19);
                end
                %pour les quad quadratique incomplet
                if(tempConnec(9) == 16)
                    matConnec(tempConnec(11),:) = tempConnec(12:27);
                end
                %pour les tetra lin
                if(tempConnec(9) == 4)
                    matConnec(tempConnec(11),:) = tempConnec(12:14);
                end
                %pour les tetra quadratique
                if(tempConnec(9) == 10)
                    matConnec(tempConnec(11),:) = tempConnec(12:21);
                end
            end
            tline = fgetl(fid);
        end
        tline = fgetl(fid);
    end
    %%dirichlet x
    if(size(findstr(tline, 'Frictionless Supports X')) ~= [0,0])
        %Sautons l'entete
        tline = fgetl(fid);tline = fgetl(fid);
        %lisons la premiere ligne
        tline = fgetl(fid);
        temp = 0;
        while ( size(findstr(tline, 'cmsel')) == [0,0])
            tempListNodes = sscanf(tline,'%f');
            for i =1: size(tempListNodes,1)
                DirichletX(temp+i,1:2) = [tempListNodes(i) , 0];% valeur de zero pour un encastrement ou symétrie
            end
            temp = size(DirichletX,1);
            tline = fgetl(fid);
        end
        tline = fgetl(fid);
    end
    %%dirichlet y
    if(size(findstr(tline, 'Frictionless Supports Y')) ~= [0,0])
        %Sautons l'entete
        tline = fgetl(fid);tline = fgetl(fid);
        %lisons la premiere ligne
        tline = fgetl(fid);
        temp = 0;
        while ( size(findstr(tline, 'cmsel')) == [0,0])
            tempListNodes = sscanf(tline,'%f');
            for i =1: size(tempListNodes,1)
                DirichletY(i+temp,1:2) = [tempListNodes(i) , 0];% valeur de zero pour un encastrement ou symétrie
            end
            temp = size(DirichletY,1);
            tline = fgetl(fid);
        end
    end
    %%dirichlet z
    if(size(findstr(tline, 'Rotational Supports Z')) ~= [0,0])
        %Sautons l'entete
        tline = fgetl(fid);tline = fgetl(fid);
        %lisons la premiere ligne
        tline = fgetl(fid);
        temp = 0;
        while ( size(findstr(tline, 'cmsel')) == [0,0])
            tempListNodes = sscanf(tline,'%f');
            for i =1: size(tempListNodes,1)
                DirichletZ(i+temp,1:2) = [tempListNodes(i) , 0];% valeur de zero pour un encastrement ou symétrie
            end
            temp = temp + size(DirichletZ,1);
            tline = fgetl(fid);
        end
    end
    if(size(findstr(tline, 'Define Temperature Constraint')) ~= [0,0])
        %Sautons l'entete
        tline = fgetl(fid);
        tline = fgetl(fid);
        %lisons la premiere ligne
        tline = fgetl(fid);
        temp = 0;
        while ( size(findstr(tline, 'fcum')) == [0,0])
            tempConnecC = sscanf(tline,'%f');
            for i =1: size(tempConnecC,1)
                T_Const(temp+i,1:2) = [tempConnecC(i) , 100+273.15];% valeur de temperature en K
            end
            temp = size(T_Const,1);
            tline = fgetl(fid);
        end
        %on saute 6 lignes
        %         for i =1:6
        %             tline = fgetl(fid);
        %         end
        %         %on prend la temperatur en Celsius
        %         fprintf('La temperature imposé est en celsius, seulement une peux etre traité.\n')
        %         posDebut = findstr(tline,'=');
        %         tImp = sscanf(tline(posDebut+1:size(tline,2)-1),'%f');
        %         T_Const(:,2) = tImp;%
        %         % fin du transfert de la temperature imposé
        %         tline = fgetl(fid);
    end
    %contour, a la suite des élément 2D dans le vecteur connec
    if (size(findstr(tline, 'Define Pressure Using Surface Effect Elements')) ~= [0,0])
        %Sautons l'entete (4 lignes)
        tline = fgetl(fid);tline = fgetl(fid);tline = fgetl(fid);tline = fgetl(fid);
        %lisons la premiere ligne
        tline = fgetl(fid);
        while ( size(tline,2) ~= 2 || size(findstr(tline, '-1'),1) == 0)
            tempConnecC = sscanf(tline,'%f');
            if(dim == 2)%elements de contour en 2 dimensions
                %pour ligne quad
                if(size(tempConnecC,1) == 8)
                    matConnec(tempConnecC(1),1:3) = tempConnecC(6:8);
                end
                %pour ligne interpol lineaire
                if(size(tempConnecC,1) == 7)
                    matConnec(tempConnecC(1),1:2) = tempConnecC(6:7);
                end
            elseif(dim == 3)%elements de contour en 3 dimensions
                %pour les quad lin en
                if(size(tempConnecC,1) == 8)
                    matConnec(tempConnecC(1),:) = tempConnecC(6:13);
                end
                %pour les quad quadratique incomplet
                if(size(tempConnecC,1) == 12)
                    matConnec(tempConnecC(1),:) = tempConnecC(6:17);
                end
                %pour les triangle lin
                if(size(tempConnecC,1) == 7)
                    matConnec(tempConnecC(1),:) = tempConnecC(6:12);
                end
                %pour les triangle quadratique
                if(size(tempConnecC,1) == 10)
                    matConnec(tempConnecC(1),:) = tempConnecC(6:15);
                end
            end
            tline = fgetl(fid);
        end
    end
    if (size(findstr(tline, 'Create "Convection"')) ~= [0,0])
        %Sautons l'entete
        for i =1:4
            tline = fgetl(fid);
        end
        %lisons la premiere ligne
        tline = fgetl(fid);
        while ( size(tline,2) ~= 2 || size(findstr(tline, '-1'),1) == 0)
            tempConnecC = sscanf(tline,'%f');
            if(dim == 2)%elements de contour en 2 dimensions
                %pour ligne quad
                if(size(tempConnecC,1) == 8)
                    matConnec(tempConnecC(1),1:3) = tempConnecC(6:8);
                end
                %pour ligne interpol lineaire
                if(size(tempConnecC,1) == 7)
                    matConnec(tempConnecC(1),1:2) = tempConnecC(6:7);
                end
            elseif(dim == 3)%elements de contour en 3 dimensions
                %pour les quad lin en
                if(size(tempConnecC,1) == 8)
                    matConnec(tempConnecC(1),:) = tempConnecC(6:13);
                end
                %pour les quad quadratique incomplet
                if(size(tempConnecC,1) == 12)
                    matConnec(tempConnecC(1),:) = tempConnecC(6:17);
                end
                %pour les triangle lin
                if(size(tempConnecC,1) == 7)
                    matConnec(tempConnecC(1),:) = tempConnecC(6:12);
                end
                %pour les triangle quadratique
                if(size(tempConnecC,1) == 10)
                    matConnec(tempConnecC(1),:) = tempConnecC(6:15);
                end
            end
            tline = fgetl(fid);
        end
        tline = fgetl(fid);
    end
    tline = fgetl(fid);
end
%fin du fichier
fclose(fid);
if (size(findstr(typeAnal, 'therm')) ~= [0,0] )
    DirichletX = T_Const;
    DirichletY = [];
    DirichletZ = [];
end
%if ( dim == 2)
%    DirichletZ = [];
%end
%pause