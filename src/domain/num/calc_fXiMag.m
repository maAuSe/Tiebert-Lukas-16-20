function fXiMag = calc_fXiMag(nF,fXi)


   fXiMag = zeros(1,nF);
   
   for jF = 1:nF
    fXiMag(jF) = sqrt(fXi(1,jF).^2 + fXi(2,jF).^2);            
   end

   
end



