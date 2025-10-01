function data = get(This,c1,c2)


   if ~haselements(This)
      data = [];
      return
   end

   if nargin == 1
      c1 = [];
      c2 = [];
   end

   if nargin == 2
      c2 = [];
   end

   % return required components
   if isempty(c1) && isempty(c2) % return all data
      data = This.data;
      return
   end

   if This.isscalar
      error('A scalar does not have components to return.')   
   end

   if isempty(c2)
      if This.isvector
         tmp = This.data;
         data = tmp(c1,:);
         return
      else
         % return tensor slice
         implementationnotfinished();
         tmp = reshape(This.data,This.dim,[]);
         data = tmp(c1,:)';
         return
      end
   end

   if isempty(c1)
      if This.isvector
         error('A vector does not have a second component to return.')
      else
         % return tensor slice
         implementationnotfinished();
         tmp = reshape(This.data,This.dim,[]);
         data = tmp(:,c2);
         return
      end
   end

   if This.isvector
      error('A vector is not indexed by two components.')
   else
      implementationnotfinished();
      tmp = reshape(This.data,This.dim,[]);
      data = tmp(c1,c2);
      return
   end


end