function c = mult_t_s(a,na,b,nb,dim)
% a is a tensor array of physical dimension dim
% na is the number of elements in na
% b is a scalar array 
% nb is the number of elements in nb
% PRE: na == nb || na==1 || nb==1

   els = eldsize(2,dim);

   if na==1 && nb==1
      if dim == 1
         c = zeros(els,1);
         c(1) = a(1) * b;
         return
      elseif dim == 2
         c = zeros(els,1);
         c(1) = a(1) * b;
         c(2) = a(2) * b;
         c(3) = a(3) * b;
         c(4) = a(4) * b;
         return
      else
         c = zeros(els,1);
         c(1) = a(1) * b;
         c(2) = a(2) * b;
         c(3) = a(3) * b;
         c(4) = a(4) * b;
         c(5) = a(5) * b;
         c(6) = a(6) * b;
         c(7) = a(7) * b;
         c(8) = a(8) * b;
         c(9) = a(9) * b;
         return
      end   
   end
   
   
   
   if nb==1
      if dim == 1
         c = zeros(els,na);
         for ii=1:na
            c(1,ii) = a(1,ii) * b;
         end
         return
      elseif dim == 2
         c = zeros(els,na);
         for ii=1:na
            c(1,ii) = a(1,ii) * b;
            c(2,ii) = a(2,ii) * b;
            c(3,ii) = a(3,ii) * b;
            c(4,ii) = a(4,ii) * b;            
         end
         return
      else
         c = zeros(els,na);
         for ii=1:na
            c(1,ii) = a(1,ii) * b;
            c(2,ii) = a(2,ii) * b;
            c(3,ii) = a(3,ii) * b;
            c(4,ii) = a(4,ii) * b;
            c(5,ii) = a(5,ii) * b;
            c(6,ii) = a(6,ii) * b;
            c(7,ii) = a(7,ii) * b;
            c(8,ii) = a(8,ii) * b;
            c(9,ii) = a(9,ii) * b;
         end
         return
      end
   end
   
   
   
   if na==1
      if dim == 1
         c = zeros(els,nb);
         for ii=1:nb
            c(1,ii) = a(1) * b(ii);
         end
         return
      elseif dim == 2
         c = zeros(els,nb);
         for ii=1:nb
            c(1,ii) = a(1) * b(ii);
            c(2,ii) = a(2) * b(ii);
            c(3,ii) = a(3) * b(ii);
            c(4,ii) = a(4) * b(ii);            
         end
         return
      else
         c = zeros(els,nb);
         for ii=1:nb
            c(1,ii) = a(1) * b(ii);
            c(2,ii) = a(2) * b(ii);
            c(3,ii) = a(3) * b(ii);
            c(4,ii) = a(4) * b(ii);
            c(5,ii) = a(5) * b(ii);
            c(6,ii) = a(6) * b(ii);
            c(7,ii) = a(7) * b(ii);
            c(8,ii) = a(8) * b(ii);
            c(9,ii) = a(9) * b(ii);            
         end
         return
      end   
   end
   
   
   
   if na == nb
      if dim == 1
         c = zeros(els,nb);
         for ii=1:na
            c(1,ii) = a(1,ii) * b(ii);
         end
         return
      elseif dim == 2
         c = zeros(els,nb);
         for ii=1:na
            c(1,ii) = a(1,ii) * b(ii);
            c(2,ii) = a(2,ii) * b(ii);
            c(3,ii) = a(3,ii) * b(ii);
            c(4,ii) = a(4,ii) * b(ii);            
         end
         return
      else
         c = zeros(els,nb);
         for ii=1:na
            c(1,ii) = a(1,ii) * b(ii);
            c(2,ii) = a(2,ii) * b(ii);
            c(3,ii) = a(3,ii) * b(ii);
            c(4,ii) = a(4,ii) * b(ii);
            c(5,ii) = a(5,ii) * b(ii);
            c(6,ii) = a(6,ii) * b(ii);
            c(7,ii) = a(7,ii) * b(ii);
            c(8,ii) = a(8,ii) * b(ii);
            c(9,ii) = a(9,ii) * b(ii);            
         end
         return
      end
   end
   
   
   c = [];
   

end