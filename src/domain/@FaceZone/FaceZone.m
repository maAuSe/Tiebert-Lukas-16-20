classdef FaceZone < RangeZone



methods (Access=public)

    
	function This = FaceZone(zoneID,dom,range)
        
   	TypeChecker.require_validrange(range,1,dom.nF);
      This = This@RangeZone(zoneID,dom,range);

   end % FaceZone CONSTRUCTOR            
   
   cdata = getcoords(This)       % CONCRETIZED
   cdata = getallcoords(This)    % CONCRETIZED
    
end % METHODS




end









