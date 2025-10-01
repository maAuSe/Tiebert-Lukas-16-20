function [datac rank dim] = pd_dot(a,b)
% PRE: isa(a,'PhysicalData') and isa(b,'PhysicalData')

   rank = a.rank + b.rank - 2;

   if a.isscalar %---------------------------------------------------------
      % [s • s]; [s • v] ; [s • t];
      error('The dot product of a scalar with another quantity is not supported.')      
   end
   
   if a.isvector %---------------------------------------------------------
      % [v • s]; v • v ; v • t;
      dim = a.dim;
      if b.isscalar
         error('The dot product of a vector with a scalar is not supported.')      
      elseif b.isvector
         if dim ~= b.dim
            error('Dimensions must be equal for vector vector dot product');
         end         
         datac = dot_v_v(get(a),getelcount(a),get(b),getelcount(b),dim);
      else % tensor
         if dim ~= b.dim
            error('Dimensions must be equal for vector tensor dot product');
         end         
         datac = dot_v_t(get(a),getelcount(a),get(b),getelcount(b),dim);
      end
      if isempty(datac)
         error('Attempt to dot product data of incompatible size.')
      end            
      return
   end
   
   if a.istensor %---------------------------------------------------------
      % [t • s]; t • v ; t • t;
      dim = a.dim;
      if b.isscalar
         error('The dot product of a tensor with a scalar is not supported.')      
      elseif b.isvector
         if dim ~= b.dim
            error('Dimensions must be equal for tensor vector dot product');
         end         
         datac = dot_t_v(get(a),getelcount(a),get(b),getelcount(b),dim);
      else % tensor
         if dim ~= b.dim
            error('Dimensions must be equal for tensor tensor dot product');
         end         
         datac = dot_t_t(get(a),getelcount(a),get(b),getelcount(b),dim);
      end
      if isempty(datac)
         error('Attempt to dot product data of incompatible size.')
      end        
      return
   end

end




