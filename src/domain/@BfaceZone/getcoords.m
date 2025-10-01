function cdata = getcoords(This)

   dom = This.dom;
         
   fCoord = dom.fCoord;          
   
   % Always returns new array
   range = This.range;
   cdata = fCoord(:,range(1):range(2));

end