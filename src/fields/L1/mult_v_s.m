function c = mult_v_s(a,na,b,nb,dim)
% a is a vector array of physical dimension dim
% na is the number of elements in na
% b is a scalar array 
% nb is the number of elements in nb
% PRE: na == nb || na==1 || nb==1

   els = eldsize(1,dim);

   if na==1 && nb==1
      if dim == 1
         c = zeros(els,1);
         c(1) = a(1) * b;                  
         return
      elseif dim == 2
         c = zeros(els,1);
         c(1) = a(1) * b;
         c(2) = a(2) * b;         
         return
      else
         c = zeros(els,1);
         c(1) = a(1) * b;
         c(2) = a(2) * b;
         c(3) = a(3) * b;
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
         end
         return
      else
         c = zeros(els,na);
         for ii=1:na
            c(1,ii) = a(1,ii) * b;
            c(2,ii) = a(2,ii) * b;
            c(3,ii) = a(3,ii) * b;
         end
         return
      end
   end
   
   
   
   if na==1
      if dim == 1
         c = zeros(int32(1),nb);
         for ii=1:nb
            c(1,ii) = a(1) * b(ii);         
         end
         return
      elseif dim == 2
         c = zeros(int32(2),nb);
         for ii=1:nb
            c(1,ii) = a(1) * b(ii);
            c(2,ii) = a(2) * b(ii);
         end
         return
      else
         c = zeros(int32(3),nb);
         for ii=1:nb
            c(1,ii) = a(1) * b(ii);
            c(2,ii) = a(2) * b(ii);
            c(3,ii) = a(3) * b(ii);
         end
         return
      end   
   end
   
   
   
   if na == nb
      if dim == 1
         c = zeros(els,nb);
         for ii=1:nb
            c(1,ii) = a(1,ii) * b(ii);         
         end
         return
      elseif dim == 2
         c = zeros(els,nb);
         for ii=1:nb
            c(1,ii) = a(1,ii) * b(ii);
            c(2,ii) = a(2,ii) * b(ii);
         end
         return
      else
         c = zeros(els,nb);
         for ii=1:nb
            c(1,ii) = a(1,ii) * b(ii);
            c(2,ii) = a(2,ii) * b(ii);
            c(3,ii) = a(3,ii) * b(ii);
         end
         return
      end
   end
   
   
   c = [];
   

end