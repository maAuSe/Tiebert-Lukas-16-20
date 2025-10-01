function display(This)
    
   if This.isscalar
      fprintf('\n Scalar Value:')
      if This.getelcount() == 1
         fprintf(['\n  ' '% 22.15e  '],This.data);
      else
         fprintf(['\n  ' '%g elements present '],This.getelcount());
      end
   elseif This.isvector
      fprintf('\n %gD Vector Value:',This.dim)
      if This.getelcount() == 1
         tmp = reshape(This.data,This.dim,[]);
         for ii = 1:This.dim
            fprintf(['\n  ' '% 22.15e  '],tmp(ii));
         end
      else
         fprintf(['\n  ' '%g elements present '],This.getelcount());
      end
   else % istensor
      fprintf('\n %gD Tensor Value:',This.dim)
      if This.getelcount() == 1
         tmp = reshape(This.data,This.dim,This.dim);
         for ii = 1:This.dim
            fprintf(['\n  ' repmat('% 22.15e   ',1,This.dim)],tmp(ii,:));
         end
      else
         fprintf(['\n  ' '%g elements present '],This.getelcount());
      end
   end
   fprintf('\n\n');
        
end