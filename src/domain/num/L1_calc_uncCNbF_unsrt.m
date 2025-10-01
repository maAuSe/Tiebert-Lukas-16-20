function [cNbFLoc uncCNbF nConnections] = L1_calc_uncCNbF_unsrt(maxnoCNbF,...
   nF,fNbCLoc,fNbC,nC)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% FUNCTION [cNbFLoc uncCNbF pCconn] = L1_calc_puncCNbF_unsrt(nPc,nF,fNbCLoc,fNbC)
%
% Calculates the topology information for indexing from physical cells to
% faces on an unstructured mesh. Input paremeters are:
%  nPc     : number of cells
%  nF      : number of faces
%  fNbCLoc      : fixed number of neighbors
%  fNbC         : contains the cell neighbors of all faces
%
% NOTE: Keep in mind the ordering of the data. Faces are not ordered cockwise or
%       counterclockwise for each cell. Another routine takes care of this!!!
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


cNbFLoc = int32([1 + maxnoCNbF*( (1:1:nC+1) - 1)]);
uncCNbF = int32(zeros(1,maxnoCNbF*nC));


for jF = 1:nF
   c1 = fNbC(fNbCLoc*(jF-1) + 1);
   basePtr = maxnoCNbF*(c1-1) + 1;
   while uncCNbF(basePtr) ~= 0
      basePtr = basePtr + 1;
   end
   uncCNbF(basePtr) = jF;

   c2 = fNbC(fNbCLoc*(jF-1) + 2);
   basePtr = maxnoCNbF*(c2-1) + 1;
   while uncCNbF(basePtr) ~= 0
      basePtr = basePtr + 1;
   end
   uncCNbF(basePtr) = jF;

end

nConnections = 2*nF;
      

end

