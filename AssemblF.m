function FG=AssemblF(FG,fe,index)
%----------------------------------------------------------
%  Purpose:
%     Assembly of element matrices into the system matrix &
%     Assembly of element vectors into the system vector
%
%  Variable Description:
%     KG - system matrix
%     FG - system vector
%     ke  - element matrix
%     fe  - element vector
%     index - d.o.f. vector associated with an element
%-----------------------------------------------------------
edof = length(index);

for i=1:edof
    ii=index(i);
    FG(ii)=FG(ii) +fe(i);
end
    