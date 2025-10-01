function c = mult_s_v(a,na,b,nb,dim)
% a is a scalar array
% na is the number of elements in na
% b is a vector array of physical dimension dim
% nb is the number of elements in nb
% PRE: na == nb || na==1 || nb==1

   els = eldsize(1,dim);

   if na==1 && nb==1
      if dim == 1
         c = zeros(els,1);
         c(1) = a * b(1);                  
         return
      elseif dim == 2
         c = zeros(els,1);
         c(1) = a * b(1);
         c(2) = a * b(2);         
         return
      else
         c = zeros(els,1);
         c(1) = a * b(1);
         c(2) = a * b(2);
         c(3) = a * b(3);
         return
      end   
   end
   
   
   
   if nb==1
      if dim == 1
         c = zeros(els,na);
         for ii=1:na
            c(1,ii) = a(ii) * b(1);         
         end
         return
      elseif dim == 2
         c = zeros(els,na);
         for ii=1:na
            c(1,ii) = a(ii) * b(1);
            c(2,ii) = a(ii) * b(2);
         end
         return
      else
         c = zeros(els,na);
         for ii=1:na
            c(1,ii) = a(ii) * b(1);
            c(2,ii) = a(ii) * b(2);
            c(3,ii) = a(ii) * b(3);
         end
         return
      end
   end
   
   
   
   if na==1
      if dim == 1
         c = zeros(els,nb);
         for ii=1:nb
            c(1,ii) = a * b(1,ii);         
         end
         return
      elseif dim == 2
         c = zeros(els,nb);
         for ii=1:nb
            c(1,ii) = a * b(1,ii);
            c(2,ii) = a * b(2,ii);
         end
         return
      else
         c = zeros(els,nb);
         for ii=1:nb
            c(1,ii) = a * b(1,ii);
            c(2,ii) = a * b(2,ii);
            c(3,ii) = a * b(3,ii);
         end
         return
      end   
   end
   
   
   
   if na == nb
      if dim == 1
         c = zeros(els,nb);
         for ii=1:na
            c(1,ii) = a(ii) * b(1,ii);         
         end
         return
      elseif dim == 2
         c = zeros(els,nb);
         for ii=1:na
            c(1,ii) = a(ii) * b(1,ii);
            c(2,ii) = a(ii) * b(2,ii);
         end
         return
      else
         c = zeros(els,nb);
         for ii=1:na
            c(1,ii) = a(ii) * b(1,ii);
            c(2,ii) = a(ii) * b(2,ii);
            c(3,ii) = a(ii) * b(3,ii);
         end
         return
      end
   end
   
   
   c = [];
   

end