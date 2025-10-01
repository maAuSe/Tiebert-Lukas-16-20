function [cVol cCoord] = calc_gcVol_gcCoord_mirror(nPc,nC,...
   cVol,cCoord,fCoord,...
   fNbCLoc,fNbC,cNbFLoc,cNbF)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Calculates the volumina and cell centroid for mirror ghost cells.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fNbC = reshape(fNbC,fNbCLoc,[]);

for jC = nPc+1:nC
   jCNbFLoc = cNbFLoc(jC);
   jF = cNbF(jCNbFLoc);
   jFNbC1 = fNbC(1,jF);
   cCoord(1,jC) = fCoord(1,jF) + (fCoord(1,jF) - cCoord(1,jFNbC1));
   cCoord(2,jC) = fCoord(2,jF) + (fCoord(2,jF) - cCoord(2,jFNbC1));      
   cVol(jC) = cVol(jFNbC1);
end


end



