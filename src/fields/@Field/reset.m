function reset(This,eldata,c1,c2)
% resets all elements in the value object with an element eldata that is
% copied along the whole current number of elements in the value object.


TypeChecker.requireclass(eldata,'double');

% Set empty value
if isempty(eldata)
   This.data = eldata;
   This.elcount = int32(0);
   return
end

if nargin == 2
   c1 = [];
   c2 = [];
end

if nargin == 3
   c2 = [];
end

if isempty(c1) && isempty(c2)
   if int32(size(eldata,1)) ~= This.elsize
      error('The supplied data does not match the storage requirements for scalars/vector/tensors.')
   end   
   This.data = repmat(eldata,1,This.elcountzone);
   This.elcount = This.elcountzone;   
   return
end

implementationnotfinished();


end