function fXiLambda = calc_fXiLambda(...
      nF,...
      fXi,fTangent,cCoord,vCoord,...
      fNbCLoc,fNbC,fNbVLoc,fNbV)
% Lambda is fraction of the node connector line starting from N1 so that
% the intersection point with the face is obtained.


   fNbC = reshape(fNbC,fNbCLoc,[]);
   fNbV = reshape(fNbV,fNbVLoc,[]);
   
   fXiLambda = double(zeros(1,nF));
   
	for jF=1:nF
    
      nbC1 = fNbC(1,jF);      
      nbV1 = fNbV(1,jF);
      
      fXiLambda(jF) = (...
         fTangent(2,jF)*(cCoord(1,nbC1) - vCoord(1,nbV1)) - ...
         fTangent(1,jF)*(cCoord(2,nbC1) - vCoord(2,nbV1))   ...
         )/(...
         fXi(2,jF)*fTangent(1,jF) - fXi(1,jF)*fTangent(2,jF) ...
         );
       
   end
   
   
   
   
   
end



