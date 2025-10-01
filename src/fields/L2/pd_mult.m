function [datac rank dim] = pd_mult(a,b)
% PRE: isa(a,'PhysicalData') and isa(b,'PhysicalData')

   rank = a.rank + b.rank;

   if a.isscalar %---------------------------------------------------------
      % s s ; s v; s t
      dim = b.dim;
      if b.isscalar         
         datac = mult_s_s(get(a),getelcount(a),get(b),getelcount(b));
      elseif b.isvector         
         datac = mult_s_v(get(a),getelcount(a),get(b),getelcount(b),dim);
      else % tensor         
         datac = mult_s_t(get(a),getelcount(a),get(b),getelcount(b),dim);
      end
      if isempty(datac)
         error('Attempt to multiply data of incompatible size.')
      end      
      return
   end
   
   if a.isvector %---------------------------------------------------------
      % v s ; v v; [v t]
      dim = a.dim;
      if b.isscalar         
         datac = mult_v_s(get(a),getelcount(a),get(b),getelcount(b),dim);
      elseif b.isvector
         if dim ~= b.dim
            error('Dimensions must be equal for vector vector dyadic product');
         end         
         datac = mult_v_v(get(a),getelcount(a),get(b),getelcount(b),dim);
      else % tensor
         error('The vector tensor product is not supported')
      end
      if isempty(datac)
         error('Attempt to multiply data of incompatible size.')
      end      
      return
   end
   
   if a.istensor %---------------------------------------------------------
      % t s ; [t v]; [t t]
      dim = a.dim;
      if b.isscalar         
         datac = mult_t_s(get(a),getelcount(a),get(b),getelcount(b),dim);
      elseif b.isvector
         error('The tensor vector product is not supported')
      else % tensor
         error('The tensor tensor product is not supported')
      end
      if isempty(datac)
         error('Attempt to multiply data of incompatible size.')
      end      
      return      
   end

end