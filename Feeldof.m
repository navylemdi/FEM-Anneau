function index=Feeldof(nd,nnel,ddln)
%----------------------------------------------------------
%  Purpose:
%     Compute system dofs associated with each element 
%
%  Synopsis:
%     [index]=feeldof(nd,nnel,ddln)
%
%  Variable Description:
%     index - system dof vector associated with element "iel"
%     nd - element number whose system dofs are to be determined
%     nnel - number of nodes per element
%     ddln - number of dofs per node 
%-----------------------------------------------------------
 
 edof = nnel*ddln;
   k=0;
   for i=1:nnel
     start = (nd(i)-1)*ddln;
       for j=1:ddln;
         k=k+1;
         index(k)=start+j;
       end
   end

 
