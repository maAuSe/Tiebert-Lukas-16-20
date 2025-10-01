function c = dot_v_v(a,na,b,nb,dim)
% a is a vector array of physical dimension dim
% na is the number of elements in na
% b is a vector array of physical dimension dim
% nb is the number of elements in nb
% PRE: na == nb || na==1 || nb==1

   if na==1 && nb==1
      if dim == 1
         c = zeros(1,1);
         c = a(1)*b(1);
         return
      elseif dim == 2
         c = zeros(1,1);
         c = a(1)*b(1) + a(2)*b(2);
         return
      else
         c = zeros(1,1);
         c = a(1)*b(1) + a(2)*b(2) + a(3)*b(3);
         return
      end   
   end
   
   
   
   if nb==1
      if dim == 1
         c = zeros(1,na);
         for ii=1:na
            c(ii) = a(1,ii)*b(1);
         end
         return
      elseif dim == 2
         c = zeros(1,na);
         for ii=1:na
            c(ii) = a(1,ii)*b(1) + a(2,ii)*b(2);
         end
         return
      else
         c = zeros(1,na);
         for ii=1:na
            c(ii) = a(1,ii)*b(1) + a(2,ii)*b(2) + a(3,ii)*b(3);
         end
         return
      end
   end
   
   
   
   if na==1
      if dim == 1
         c = zeros(1,nb);
         for ii=1:na
            c(ii) = a(1)*b(1,ii);
         end
         return
      elseif dim == 2
         c = zeros(1,nb);
         for ii=1:na
            c(ii) = a(1)*b(1,ii) + a(2)*b(2,ii);
         end
         return
      else
         c = zeros(1,nb);
         for ii=1:na
            c(ii) = a(1)*b(1,ii) + a(2)*b(2,ii) + a(3)*b(3,ii);
         end
         return
      end
   end
   
   
   
   if na==nb
      if dim == 1
         c = zeros(1,nb);
         for ii=1:nb
            c(ii) = a(1,ii)*b(1,ii);
         end
         return
      elseif dim == 2
         c = zeros(1,nb);
         for ii=1:nb
            c(ii) = a(1,ii)*b(1,ii) + a(2,ii)*b(2,ii);
         end
         return
      else
         c = zeros(1,nb);
         for ii=1:nb
            c(ii) = a(1,ii)*b(1,ii) + a(2,ii)*b(2,ii) + a(3,ii)*b(3,ii);
         end
         return
      end
   end
   
   
   c = [];
   

end