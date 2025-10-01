classdef BfaceZone < FaceZone


   
methods (Access=public)

    
	function This = BfaceZone(zoneID,dom,range)
        
      TypeChecker.require_validrange(range,dom.nIf+1,dom.nF);
      This = This@FaceZone(zoneID,dom,range);

   end % BfaceZone CONSTRUCTOR            

   cdata = getcoords(This)       % OVERLOADED
   cdata = getallcoords(This)    % OVERLOADED
        
end % METHODS




end









