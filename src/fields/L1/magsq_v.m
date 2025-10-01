function b = magsq_v(a,na,dim)

   els = eldsize(0,dim);

   if na==1
      if dim == 1
         b = zeros(els,1);
         b(1) = a*a;
         return
      elseif dim == 2
         b = zeros(els,1);         
         b(1) = a(1)*a(1) + a(2)*a(2);
      else
         b = zeros(els,1);         
         b(1) = a(1)*a(1) + a(2)*a(2) + a(3)*a(3);
      end
      return
   else
      if dim == 1
         b = zeros(els,na);
         for ii=1:na
            b(ii) = a(ii)*a(ii);
         end
      elseif dim == 2
         b = zeros(els,na);         
         for ii=1:na
            b(ii) = a(1,ii)*a(1,ii) + a(2,ii)*a(2,ii);
         end            
      else
         b = zeros(els,na);         
         for ii=1:na
            b(ii) = a(1,ii)*a(1,ii) + a(2,ii)*a(2,ii) + a(3,ii)*a(3,ii);
         end
      end
      return      
   end
   
end