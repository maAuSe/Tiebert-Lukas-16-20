function [datac rank dim] = pd_subtr(a,b)
% PRE: isa(a,'PhysicalData') and isa(b,'PhysicalData')

   rank = a.rank;

   if a.isscalar %---------------------------------------------------------
      % s - s     
      dim = b.dim;
      if b.isscalar
         datac = subtr_s_s(get(a),getelcount(a),get(b),getelcount(b));         
      elseif b.isvector         
         error('Cannot subtract a scalar and a vector.')
      else % tensor
         error('Cannot subtract a scalar and a tensor.')
      end
      if isempty(datac)
         error('Attempt to subtract data of incompatible size.')
      end
      return
   end
   
   if a.isvector %---------------------------------------------------------
      % v - v      
      dim = a.dim;
      if b.isscalar
         error('Cannot subtract a vector and a scalar.')
      elseif b.isvector
         if dim ~= b.dim
            error('Dimensions must be equal when subtracting two vectors.');
         end
         datac = subtr_v_v(get(a),getelcount(a),get(b),getelcount(b),dim);
      else % tensor
         error('Cannot subtract a vector and a tensor.')
      end
      if isempty(datac)
         error('Attempt to subtract data of incompatible size.')
      end
      return
   end
   
   if a.istensor %---------------------------------------------------------
      % t - t      
      dim = a.dim;
      if b.isscalar
         error('Cannot subtract a tensor and a scalar.')
      elseif b.isvector
         error('Cannot subtract a tensor and a vector.')
      else % tensor
         if dim ~= b.dim
            error('Dimensions must be equal when subtracting two tensors.');
         end
         datac = subtr_t_t(get(a),getelcount(a),get(b),getelcount(b),dim);
      end
      if isempty(datac)
         error('Attempt to subtract data of incompatible size.')
      end
      return      
   end
end






