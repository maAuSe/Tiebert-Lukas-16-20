function cdata = getcoords(This)

   dom = This.dom;
      
   nF     = dom.nF;
   fCoord = dom.fCoord;
    
   elc = getelcount(This);   
      
   if elc == nF
      % Return whole array
      cdata = fCoord;
   else
      % Return part of the array
      range = This.range;
      cdata = fCoord(:,range(1):range(2));
   end


end