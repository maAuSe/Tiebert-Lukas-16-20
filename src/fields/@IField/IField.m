classdef IField < handle
% Interface for fields

methods (Access=public,Abstract=true)
   
   
   DomObj = getdomain(This)      % returns a reference to the domain
   ZoneObj = getzone(This)       % returns a reference to the zone
   coordsdata = getcoords(This)  % returns the point coordinates of the zone 
   n = getelcountzone(This)      % returns the number of points in the zone
   FieldObj = restrictto(This,zoneObj) % returns new field on subzone
    
   
end % ABSTRACT METHODS


end