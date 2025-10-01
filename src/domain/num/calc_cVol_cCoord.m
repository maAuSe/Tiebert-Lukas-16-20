function [cVol cCoord] = calc_cVol_cCoord(nC,nPc,...
   vCoord,fCoord,...
   cNbFLoc,cNbF,fNbVLoc,fNbV,fNbCLoc,fNbC)

   
cVol   = zeros(1,nC);
cCoord = zeros(2,nC);


[cVol cCoord] = calc_pCVol_pCCoord(nPc,...
   cVol,cCoord,...
   cNbFLoc,cNbF,fNbVLoc,fNbV,vCoord);


[cVol cCoord] = calc_gcVol_gcCoord_mirror(nPc,nC,...
   cVol,cCoord,fCoord,...
   fNbCLoc,fNbC,cNbFLoc,cNbF);

end



