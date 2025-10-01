function cdata = getcoords(This)

   dom = This.dom;
      
   nV     = dom.nV;
   vCoord = dom.vCoord;
    
   elc = getelcount(This);   
      
   if elc == nV 
      % Return whole array
      cdata = vCoord;
   else
      % Return part of the array
      range = This.range;
      cdata = vCoord(:,range(1):range(2));
   end


end