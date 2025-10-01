function [cVol cCoord] = calc_pCVol_pCCoord(nPc,...
   cVol,cCoord,...
   cNbFLoc,cNbF,fNbVLoc,fNbV,vCoord)



    fNbV = reshape(fNbV,fNbVLoc,[]);
    
    jFStartLoc = cNbFLoc(1);
    for jC = 1:nPc

        nextJfStartLoc = cNbFLoc(jC+1);
        
        refV = fNbV(1,cNbF(jFStartLoc));
        refVX = vCoord(1,refV);
        refVY = vCoord(2,refV);
        
        cVol_tmp = 0;
        Xmom_tmp = 0;
        Ymom_tmp = 0;

        % Loop over all connected faces, except the first
        % Note that the last face cannot be excluded because we have no
        % control wether the refV is V1 or V2 of the second face.
        for jFLoc= jFStartLoc:nextJfStartLoc-1

            jF = cNbF(jFLoc);            
            V1 = fNbV(1,jF);
            V2 = fNbV(2,jF);            
            d1X = vCoord(1,V1) - refVX;
            d1Y = vCoord(2,V1) - refVY;
            d2X = vCoord(1,V2) - refVX;
            d2Y = vCoord(2,V2) - refVY;
            % temp is the area of the triangle Vref->V1->V2
            triArea = abs(0.5*(d2X*d1Y - d2Y*d1X)); % (2009-10-19) We changed sign of this            
            triX = (d1X + d2X)/3;
            triY = (d1Y + d2Y)/3;
         
            cVol_tmp = cVol_tmp + triArea;

            Xmom_tmp = Xmom_tmp + triX*triArea;
            Ymom_tmp = Ymom_tmp + triY*triArea;

        end % Loop over all faces connected to current cell                

        cVol(jC) = cVol_tmp;

        cCoord(1,jC) = Xmom_tmp/cVol_tmp + refVX;
        cCoord(2,jC) = Ymom_tmp/cVol_tmp + refVY;

        jFStartLoc = nextJfStartLoc;

    end % loop over cells

end




