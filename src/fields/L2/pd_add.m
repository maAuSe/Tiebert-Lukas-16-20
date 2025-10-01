function [datac rank dim] = pd_add(a,b)
% PRE: isa(a,'PhysicalData') and isa(b,'PhysicalData')

   rank = a.rank;

   if a.isscalar %---------------------------------------------------------
      % s + s      
      dim = b.dim;
      if b.isscalar
         datac = add_s_s(get(a),getelcount(a),get(b),getelcount(b));
      elseif b.isvector         
         error('Cannot add a scalar and a vector.')
      else % tensor
         error('Cannot add a scalar and a tensor.')
      end
      if isempty(datac)
         error('Attempt to add data of incompatible size.')
      end
      return
   end
   
   if a.isvector %---------------------------------------------------------
      % v + v      
      dim = a.dim;
      if b.isscalar
         error('Cannot add a vector and a scalar.')
      elseif b.isvector
         if dim ~= b.dim
            error('Dimensions must be equal when adding two vectors.');
         end
         datac = add_v_v(get(a),getelcount(a),get(b),getelcount(b),dim);
      else % tensor
         error('Cannot add a vector and a tensor.')
      end
      if isempty(datac)
         error('Attempt to add data of incompatible size.')
      end
      return
   end
   
   if a.istensor %---------------------------------------------------------
      % t + t      
      dim = a.dim;
      if b.isscalar
         error('Cannot add a tensor and a scalar.')
      elseif b.isvector
         error('Cannot add a tensor and a vector.')
      else % tensor
         if dim ~= b.dim
            error('Dimensions must be equal when adding two tensors.');
         end
         datac = add_t_t(get(a),getelcount(a),get(b),getelcount(b),dim);
      end
      if isempty(datac)
         error('Attempt to add data of incompatible size.')
      end
      return      
   end

end


