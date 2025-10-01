classdef Value < PhysicalData & IValue
   
   
properties (GetAccess=public,SetAccess=protected)

   data     % double(elsize,elcount)
   % Holds the data
   elcount  % int32
   % Holds the number of data elements.

end % PROPERTIES



methods (Access=public,Static=true)
   
   
   
end % STATIC METHODS



methods (Access=public)

   function Obj = Value(rank,dim)
   %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   % Value object constructor. Creates an empty container for storing
   % values of the given rank and dimension of the physical problem.
   % POST:
   %  A valid Value object iff DATA meets all the following conditions:   
   %  (1) rank is one of the following 0,1,2
   %  (2) dim is one of the following 1,2,3
   %
   %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      Obj = Obj@PhysicalData(rank,dim);
      Obj = Obj@IValue;
      
      Obj.data = [];
      Obj.elcount = int32(0);

   end

   % Implementation of methods PhysicalData   
   elcount = getelcount(This)
   haselflag = haselements(This)
   
   % Implementation of methods in IValue
   data = get(This,c1,c2)
   data = extract(This)
   set(This,data,c1,c2)
   reset(This,eldata,c1,c2)
   
   % Implementation of methods in IOverloadableOperators
   c = plus(a,b)
   c = minus(a,b)
   c = mtimes(a,b)
   
   c = uminus(a)
   c = uplus(a)   
   c = ctranspose(a)
   
   % Other operations
   c = dot(a,b) % TODO: implement dot(a,b,c)
   
   c = magsq(a);
   c = mag(a);

	% Implementation of auxiliary methods for binary operators
   c = lhs_plus(a,b)
   c = rhs_plus(b,a)
   c = lhs_minus(a,b)
   c = rhs_minus(b,a)
   c = lhs_mtimes(a,b)
   c = rhs_mtimes(b,a)
   c = lhs_dot(a,b)
   c = rhs_dot(b,a)
      
   % OLD
   % objCopy = copy(obj)   
   % c = negadd(a,b)   
   % c = elinvmult(a,b)      

end % METHODS

   
   
   
end












