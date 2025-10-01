function a = calc_MaxNoCNbV(nPc,cNbVLoc)
% Calculates the maximum of the number of vertex neighbors of each cell in
% the mesh.

   a = 0;

   for jC = 1:nPc
      thisCNoNbV = cNbVLoc(jC+1) - cNbVLoc(jC);
      a = max(a,thisCNoNbV);
   end


end