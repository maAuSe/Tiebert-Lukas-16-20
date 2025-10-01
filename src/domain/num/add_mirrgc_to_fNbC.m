function fNbC = add_mirrgc_to_fNbC(fNbCLoc,fNbC,nIf,nBf,nPc)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% FUNCTION fNbC = add_mirrgc_to_fNbC(fNbC,fStart,nBf,nPc)
%
% Adds mirror ghost cells to fNbC as follows:
%  for boundary face jBf
%     fNbC1 remains
%     fNbC2 becomes nPc + jBf
%  end
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


   for jBf = 1:1:nBf

      jF = nIf + jBf;

      jFNbCLoc = fNbCLoc*(jF-1)+1;
      fNbC(jFNbCLoc+1) = nPc + jBf;

   end


end

