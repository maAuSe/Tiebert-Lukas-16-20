classdef VertZone < RangeZone


methods (Access=public)

    
	function This = VertZone(zoneID,dom,range)
        
   	TypeChecker.require_validrange(range,1,dom.nV);
      This = This@RangeZone(zoneID,dom,range);

	end % VertZone CONSTRUCTOR            

   cdata = getcoords(This)       % CONCRETIZED
   cdata = getallcoords(This)    % CONCRETIZED
    
end % METHODS




end









