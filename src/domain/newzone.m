function Obj = newzone(name,dom,range,kind)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This factory routine provides a wrapper for easy creation of zone objects
% defined on a domain.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   
switch kind
   case 'Vert'
      Obj = VertZone(name,dom,range);
   case 'Face'
      Obj = FaceZone(name,dom,range);
	case 'Bface'
      Obj = BfaceZone(name,dom,range);
   case 'Cell'
      Obj = CellZone(name,dom,range);
   otherwise
      error('Choose kind from {''Vert'',''Face'',''Bface'',''Cell''}');            
end
   
   
end