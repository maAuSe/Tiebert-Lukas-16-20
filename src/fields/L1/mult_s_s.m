function c = mult_s_s(a,na,b,nb)
% a is a scalar array
% na is the number of elements in na
% b is a scalar array
% nb is the number of elements in nb
% PRE: na == nb || na==1 || nb==1

   els = eldsize(0,1);
   
   if na==1 && nb==1
      c = zeros(els,1);
      c = a * b;
      return
   end
   
   
   
   if nb==1
      c = zeros(els,na);
      for ii=1:na
         c(ii) = a(ii) * b;
      end   
      return
   end
   
   
   
   if na==1
      c = zeros(els,nb);
      for ii=1:nb
         c(ii) = a * b(ii);
      end
      return
   end
   
   
   
   if na == nb
      c = zeros(els,na);
      for ii=1:na
         c(ii) = a(ii) * b(ii);
      end
      return
   end
   
   
   c = [];
   
   
end