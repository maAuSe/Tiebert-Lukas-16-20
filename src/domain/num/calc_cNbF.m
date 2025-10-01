function [cNbFLoc cNbF] = calc_cNbF(nC,nPc,nF,maxnoCNbF,fNbCLoc,fNbC,fNbVLoc,fNbV)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% FUNCTION [cNbFLoc cNbF] = calc_cNbF(...
%                 nC,nF,maxnoCNbF,fNbCLoc,fNbC,fNbVLoc,fNbV)
%
% Calculates the topology information for indexing from cells to faces
% on an unstructured mesh. Input paremeters are:
%  nC       : number of cells
%  nF       : number of faces
%  maxnoCNbF: max no face neighbors connecting to a single cell
%  fNbCLoc  : fixed number of neighbors
%  fNbC     : contains the cell neighbors of all faces
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    



%==========================================================================
% COMPUTE cNbF and SORT pCell entries appropriately
%==========================================================================
[cNbFLoc uncCNbF nConnections] = L1_calc_uncCNbF_unsrt(...
	maxnoCNbF,nF,fNbCLoc,fNbC,nC);
uncCNbF = sort_uncCNbF(nPc,cNbFLoc,uncCNbF,fNbVLoc,fNbV);




%==========================================================================
% COMPRESS uncCNbF into cNbF
%==========================================================================
cNbF = int32(zeros(1,nConnections));
[cNbFLoc cNbF] = L1_compress_ndxpair(cNbFLoc,uncCNbF,cNbF);
   

end

