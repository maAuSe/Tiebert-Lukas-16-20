classdef CellZone < RangeZone



methods (Access=public)

    
	function This = CellZone(zoneID,dom,range)
        
      TypeChecker.require_validrange(range,1,dom.nC);
      This = This@RangeZone(zoneID,dom,range);

	end % CellZone CONSTRUCTOR            

   cdata = getcoords(This)       % CONCRETIZED
   cdata = getallcoords(This)    % CONCRETIZED
    
end % METHODS




end









