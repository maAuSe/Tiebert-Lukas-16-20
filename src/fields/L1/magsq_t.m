function b = magsq_t(a,na,dim)

   els = eldsize(0,dim);

   if na==1
      if dim == 1
         b = zeros(els,1);
         b(1) = 0.5*(a*a);
         return
      elseif dim == 2
         b = zeros(els,1);         
         b(1) = 0.5*(a(1)*a(1) + a(2)*a(2) + ...
                     a(3)*a(3) + a(4)*a(4));
      else
         b = zeros(els,1);
         b(1) = 0.5*(a(1)*a(1) + a(2)*a(2) + a(3)*a(3) + ...
                     a(4)*a(4) + a(5)*a(5) + a(6)*a(6) + ...
                     a(7)*a(7) + a(8)*a(8) + a(9)*a(9));
      end
      return
   else
      if dim == 1
         b = zeros(els,na);
         for ii=1:na
            b(ii) = 0.5*(a(ii)*a(ii));
         end
      elseif dim == 2
         b = zeros(els,na);         
         for ii=1:na
            b(ii) = 0.5*(a(1,ii)*a(1,ii) + a(2,ii)*a(2,ii) + ...
                         a(3,ii)*a(3,ii) + a(4,ii)*a(4,ii));
         end            
      else
         b = zeros(els,na);         
         for ii=1:na
            b(ii) = 0.5*(a(1,ii)*a(1,ii) + a(2,ii)*a(2,ii) + a(3,ii)*a(3,ii) + ...
                         a(4,ii)*a(4,ii) + a(5,ii)*a(5,ii) + a(6,ii)*a(6,ii) + ...
                         a(7,ii)*a(7,ii) + a(8,ii)*a(8,ii) + a(9,ii)*a(9,ii));
         end
      end
      return      
   end
   
end