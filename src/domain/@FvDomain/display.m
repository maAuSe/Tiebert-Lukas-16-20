function display(This)

   fprintf('\nCLASS %s (%gD)',class(This),This.dim)    
   %fprintf('\n\t      dim : %8g',This.dim);
   fprintf('\n%8g Physical Cells',This.nPc);
   fprintf('\n%8g Ghost Cells',This.nGc);
   fprintf('\n%8g Faces (incl. boundary faces)',This.nF);
   fprintf('\n%8g Vertices (incl. boundary vertices)',This.nV);
   fprintf('\n\n');

end % FvMesh.display