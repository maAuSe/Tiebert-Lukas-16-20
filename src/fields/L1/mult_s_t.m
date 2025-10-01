function c = mult_s_t(a,na,b,nb,dim)
% a is a scalar array
% na is the number of elements in na
% b is a tensor array of physical dimension dim
% nb is the number of elements in nb
% PRE: na == nb || na==1 || nb==1

   els = eldsize(2,dim);

   if na==1 && nb==1
      if dim == 1
         c = zeros(els,1);
         c(1) = a * b(1);
         return
      elseif dim == 2
         c = zeros(els,1);
         c(1) = a * b(1);
         c(2) = a * b(2);
         c(3) = a * b(3);
         c(4) = a * b(4);
         return
      else
         c = zeros(els,1);
         c(1) = a * b(1);
         c(2) = a * b(2);
         c(3) = a * b(3);
         c(4) = a * b(4);
         c(5) = a * b(5);
         c(6) = a * b(6);
         c(7) = a * b(7);
         c(8) = a * b(8);
         c(9) = a * b(9);
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
            c(3,ii) = a(ii) * b(3);
            c(4,ii) = a(ii) * b(4);            
         end
         return
      else
         c = zeros(els,na);
         for ii=1:na
            c(1,ii) = a(ii) * b(1);
            c(2,ii) = a(ii) * b(2);
            c(3,ii) = a(ii) * b(3);
            c(4,ii) = a(ii) * b(4);
            c(5,ii) = a(ii) * b(5);
            c(6,ii) = a(ii) * b(6);
            c(7,ii) = a(ii) * b(7);
            c(8,ii) = a(ii) * b(8);
            c(9,ii) = a(ii) * b(9);            
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
            c(3,ii) = a * b(3,ii);
            c(4,ii) = a * b(4,ii);            
         end
         return
      else
         c = zeros(els,nb);
         for ii=1:nb
            c(1,ii) = a * b(1,ii);
            c(2,ii) = a * b(2,ii);
            c(3,ii) = a * b(3,ii);
            c(4,ii) = a * b(4,ii);
            c(5,ii) = a * b(5,ii);
            c(6,ii) = a * b(6,ii);
            c(7,ii) = a * b(7,ii);
            c(8,ii) = a * b(8,ii);
            c(9,ii) = a * b(9,ii);            
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
            c(3,ii) = a(ii) * b(3,ii);
            c(4,ii) = a(ii) * b(4,ii);            
         end
         return
      else
         c = zeros(els,nb);
         for ii=1:na
            c(1,ii) = a(ii) * b(1,ii);
            c(2,ii) = a(ii) * b(2,ii);
            c(3,ii) = a(ii) * b(3,ii);
            c(4,ii) = a(ii) * b(4,ii);
            c(5,ii) = a(ii) * b(5,ii);
            c(6,ii) = a(ii) * b(6,ii);
            c(7,ii) = a(ii) * b(7,ii);
            c(8,ii) = a(ii) * b(8,ii);
            c(9,ii) = a(ii) * b(9,ii);            
         end
         return
      end
   end
   
   
   c = [];
   

end