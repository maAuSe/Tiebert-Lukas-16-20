classdef IValue < handle
% Interface for values

methods (Access=public,Abstract=true)
   
   
   data = get(This,c1,c2)   % returns a copy of the object data
   data = extract(This)     % returns the object data and empties the object
   set(This,data,c1,c2)     % puts data in the object
   reset(This,datael,c1,c2) % changes the object data to the vlaue of datael
    
   
end % ABSTRACT METHODS


end