function b = magsq_s(a,na)
    
   els = eldsize(0,dim);
   
   if na==1
      b = zeros(els,1);
      b(1) = a*a;
      return
   else
      b = zeros(els,na);
      for ii=1:na
         b(ii) = a(ii)*a(ii);
      end
      return
   end
   
   

end