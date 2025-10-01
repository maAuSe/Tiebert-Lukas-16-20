classdef RangeZone < IZone & FvmLabEntity



properties (GetAccess=public,SetAccess=protected)

    dom           % FvDomain
    % The finite volume domain on which this zone is defined
    
    elcount       % int32
    % number of elements in this zone = 1 + range(2) - range(1);
    
    range         % int32(2x1)
    % First value indicates the first element in the range
    % Second value indicates the last element in the range
    
end % PROPERTIES




methods (Access=public)

    
    function This = RangeZone(zoneID,dom,range)
       % PRE: valid range spec
        
        This = This@IZone();
        This = This@FvmLabEntity(zoneID);
        
        TypeChecker.requireclass(dom,'FvDomain');        
        %TypeChecker.requirearraysize(range,2,1);
                
        This.dom = dom;
        This.range = int32(range);
        This.elcount = 1 + range(2) - range(1);
                
    end

       	
   delete(This)
   DomObj = getdomain(This)
   dim    = getdim(This)
   n      = getelcount(This)   
   range  = getrange(This)
           
end




end










