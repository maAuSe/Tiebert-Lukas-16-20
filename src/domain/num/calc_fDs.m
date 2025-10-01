function fDs = calc_fDs(nF,...
   fCoord,cCoord,...
   fNbCLoc,fNbC)
   

   fDs = double(zeros(2,2*nF));
   
   fNbC = reshape(fNbC,fNbCLoc,[]);
   
   jFLoc = int32(0);
   for jF = 1:nF
      
      % For owner side:
      jFNbC1 = fNbC(1,jF);
      jFLoc = jFLoc + 1;
      fDs(1,jFLoc) = fCoord(1,jF) - cCoord(1,jFNbC1);
      fDs(2,jFLoc) = fCoord(2,jF) - cCoord(2,jFNbC1);
            
      % For neighbor side:
      jFNbC2 = fNbC(2,jF);
      jFLoc = jFLoc + 1;
      fDs(1,jFLoc) = fCoord(1,jF) - cCoord(1,jFNbC2);
      fDs(2,jFLoc) = fCoord(2,jF) - cCoord(2,jFNbC2);
      
   end
   
   
end



