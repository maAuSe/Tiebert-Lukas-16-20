function [cNbVLoc cNbV] = calc_cNbV(nC,nPc,maxnoCNbV,cNbFLoc,cNbF,fNbVLoc,fNbV)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% FUNCTION [cNbVLoc cNbV] = calc_cNbV(...
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
% COMPUTE cNbVLoc and uncCNbV
%==========================================================================
% PRE: cNbF must be sorted such that consecutive face neighbors are
% physically connected.
cNbVLoc = int32([1 + maxnoCNbV*( (1:1:nC+1) - 1)]);
uncCNbV = int32(zeros(1,maxnoCNbV*nC));

fNbV = reshape(fNbV,fNbVLoc,[]);
nConnections = 0;
fStartLoc = cNbFLoc(1);

for jC = 1:nPc   % Loop pCells
   vLoc = maxnoCNbV*(jC-1) + 1;   
   fNextStartLoc = cNbFLoc(jC+1);
   
   % Add first vertex
   jFLoc = fStartLoc;
   jF = cNbF(jFLoc);
   v = fNbV(1,jF);
   
   jNextFace = cNbF(jFLoc+1);
   if v ~= fNbV(1,jNextFace) && v ~= fNbV(2,jNextFace)
   % If v is not shared with next face, take v = fNbV(2,jF), in order to
   % have a shared v.
      v = fNbV(2,jF);
   end   
   uncCNbV(vLoc) = v;   
   vLoc = vLoc + 1;
   vPrev = v;
   nConnections = nConnections + 1;      
   
   % Add one vertex for every every additional face (2D & valid for pCells)   
   for jFLoc = fStartLoc+1:fNextStartLoc-1
      jF = cNbF(jFLoc);
      v = fNbV(1,jF);
      if v == vPrev % sorted structure in 2D
         v = fNbV(2,jF);
      end
      uncCNbV(vLoc) = v;   
      vLoc = vLoc + 1;
      vPrev = v;
      nConnections = nConnections + 1;
   end
   
   fStartLoc = fNextStartLoc;
   
end


for jC = nPc+1:nC % loop gCells
   vLoc = maxnoCNbV*(jC-1) + 1;   
   fNextStartLoc = cNbFLoc(jC+1);
   
   jFLoc = fStartLoc;
   jF = cNbF(jFLoc);
   
   % Add first vertex   
   v = fNbV(1,jF);
   uncCNbV(vLoc) = v;
   nConnections = nConnections + 1;
   
   % Add second vertex   
   v = fNbV(2,jF);
   uncCNbV(vLoc + 1) = v;
   nConnections = nConnections + 1;
     
   fStartLoc = fNextStartLoc;
   
end




%==========================================================================
% COMPRESS uncCNbV into cNbV
%==========================================================================
cNbV = int32(zeros(1,nConnections));
[cNbVLoc cNbV] = L1_compress_ndxpair(cNbVLoc,uncCNbV,cNbV); % TODO: rename routine
   

end



