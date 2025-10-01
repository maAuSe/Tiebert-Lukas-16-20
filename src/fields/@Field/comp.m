function c = comp(This,c1,c2)

   if nargin == 1
      c1 = [];
      c2 = [];
   end

   if nargin == 2
      c2 = [];
   end

   % return required components
   if isempty(c1) && isempty(c2) % return all data
      c = Field(getzone(This),This.rank);
      set(c,get(This));
      return
   end

   if This.isscalar
      error('A scalar is not indexed by components.')   
   end

   if isempty(c2)
      if This.isvector
         tmp = This.data;
         data = tmp(c1,:);
         c = Field(getzone(This),0);
         set(c,data)
         return
      else
         % return tensor slice
         implementationnotfinished();
      end
   end

   if isempty(c1)
      if This.isvector
         error('A vector does not have a second component to return.')
      else
         % return tensor slice
         implementationnotfinished();
      end
   end

   if This.isvector
      error('A vector is not indexed by two components.')
   else
      implementationnotfinished();      
   end


end