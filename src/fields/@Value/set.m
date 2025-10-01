function set(This,data,c1,c2)

% Set empty value
if isempty(data)
   This.data = data;
   This.elcount = int32(0);
   return
end

TypeChecker.requireclass(data,'double');

if nargin == 2
   c1 = [];
   c2 = [];
end

if nargin == 3
   c2 = [];
end

if isempty(c1) && isempty(c2)
   if int32(size(data,1)) ~= This.elsize
      error('The supplied data does not match the storage requirements for scalars/vector/tensors.')
   end
   This.data = data;
   This.elcount = int32(size(data,2));
   return
end


implementationnotfinished();


if This.isscalar
	error('A scalar does not have components to set.')
end

if This.isvector
   if ~isempty(c2)
      error('A vector does not have a second index.')
   end
   % TODO: check datasize
   This.data(c1,:) = data;
   % Number of elements remains the same
   implementationnotfinished();
end

if This.istensor
   implementationnotfinished();
end


end