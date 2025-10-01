function cdata = getcoords(This)

   dom = This.dom;
      
   nC     = dom.nC;
   cCoord = dom.cCoord;
    
   elc = getelcount(This);
         
   if elc == nC
      % Return whole array
      cdata = cCoord;
   else
      % Return part of the array
      range = This.range;
      cdata = cCoord(:,range(1):range(2));
   end


end