function Obj = newvalue(data,elcount)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This factory routine provides a wrapper for easy creation of Value
% objects with a single element from given data. The supplied data is
% analysed to retrieve the rank and dimension. If elcount is supplied, that
% resulting value will have that number of elements.
% POST:
%  A valid Value object iff DATA is of class double and meets one of the
%  following conditions:
%  (1) a = size(data);
%      a(1) == 1 && a(2) == 1
%  (2) a = size(data);
%      (a(1) == (1||2||3) && a(2) == 1) || 
%      (a(2) == (1||2||3) && a(1) == 1)
%  (3) a = size(data)
%      a(1) == a(2) && a(1) == (1||2||3)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


   TypeChecker.requireclass(data,'double');
   TypeChecker.requirenonempty(data);
   
   if nargin == 1
      elcount = 1;
   end
   
   ds = size(data);
   
   if ds(1)==1 && ds(2)==1
      Obj = Value(0,1);
      Obj.set(repmat(data,1,elcount));
      % Note: dim is irrelevant for scalars thus we set it to 1.
   else      
      if ds(2) == 1 && (ds(1) == 3 || ds(1) == 2 || ds(1) == 1)
         % Vector in column notation
         Obj = Value(1,ds(1));
         Obj.set(repmat(data,1,elcount));
      elseif ds(1) == 1 && (ds(2) == 3 || ds(2) == 2 || ds(2) == 1)
         % Vector in row notation
         Obj = Value(1,ds(2));
         data = data';
         Obj.set(repmat(data,1,elcount));
      elseif ds(2) == 1 && (ds(1) == 4 || ds(1) == 9)
         % Tensor in column notation
         Obj = Value(2,sqrt(ds(1)));
         Obj.set(repmat(data,1,elcount));
      elseif ds(1) == 1 && (ds(2) == 4 || ds(2) == 9)
         % Tensor in row notation
         Obj = Value(2,sqrt(ds(2)));
         data = data';
         Obj.set(repmat(data,1,elcount));
      elseif ds(1) == ds(2) && (ds(1) == 3 || ds(1) == 2)
         % Tensor in matrix notation
         Obj = Value(2,ds(2));
         data = reshape(data,ds(1)*ds(1),[]);
         Obj.set(repmat(data,1,elcount));
      else
         error('The supplied %gx%g data cannot be fitted into a scalar/vector/tensor value',ds(1),ds(2))
      end
   end


end