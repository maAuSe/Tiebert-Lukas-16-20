function c = mult_v_v(a,na,b,nb,dim)
% a is a vector array of physical dimension dim
% na is the number of elements in a
% b is a vector array of physical dimension dim
% nb is the number of elements in b
% PRE: na == nb || na==1 || nb==1
   
   crank = 2;                 
   els = eldsize(crank,dim); % result is tensor of dimension dim

   if na==1 && nb==1
      if dim == 1
         c = zeros(els,1);
         c(1,1) = a(1,1)*b(1,1);
         return
      elseif dim == 2
         c = zeros(els,1);
         c(1,1) = a(1,1)*b(1,1);
         c(2,1) = a(2,1)*b(1,1);
         c(3,1) = a(1,1)*b(2,1);
         c(4,1) = a(2,1)*b(2,1);
         return
      else
         c = zeros(els,1);
         c(1,1) = a(1,1)*b(1,1);
         c(2,1) = a(2,1)*b(1,1);
         c(3,1) = a(3,1)*b(1,1);
         c(4,1) = a(1,1)*b(2,1);
         c(5,1) = a(2,1)*b(2,1);
         c(6,1) = a(3,1)*b(2,1);
         c(7,1) = a(1,1)*b(3,1);
         c(8,1) = a(2,1)*b(3,1);
         c(9,1) = a(3,1)*b(3,1);
         return
      end
   end
   
   
   if nb==1
      if dim == 1
         c = zeros(els,na);
         for ii=1:na
            c(1,ii) = a(1,ii)*b(1,1);
         end
         return
      elseif dim == 2
         c = zeros(els,na);
         for ii=1:na
            c(1,ii) = a(1,ii)*b(1,1);
            c(2,ii) = a(2,ii)*b(1,1);
            c(3,ii) = a(1,ii)*b(2,1);
            c(4,ii) = a(2,ii)*b(2,1);
         end
         return
      else
         c = zeros(els,na);
         for ii=1:na
            c(1,ii) = a(1,ii)*b(1,1);
            c(2,ii) = a(2,ii)*b(1,1);
            c(3,ii) = a(3,ii)*b(1,1);
            c(4,ii) = a(1,ii)*b(2,1);
            c(5,ii) = a(2,ii)*b(2,1);
            c(6,ii) = a(3,ii)*b(2,1);
            c(7,ii) = a(1,ii)*b(3,1);
            c(8,ii) = a(2,ii)*b(3,1);
            c(9,ii) = a(3,ii)*b(3,1);
         end
         return
      end      
   end
   
   
   if na==1
      if dim == 1
         c = zeros(els,nb);
         for ii=1:nb
            c(1,ii) = a(1,1)*b(1,ii);
         end
         return
      elseif dim == 2
         c = zeros(els,nb);
         for ii=1:nb
            c(1,ii) = a(1,1)*b(1,ii);
            c(2,ii) = a(2,1)*b(1,ii);
            c(3,ii) = a(1,1)*b(2,ii);
            c(4,ii) = a(2,1)*b(2,ii);
         end
         return
      else
         c = zeros(els,nb);
         for ii=1:nb
            c(1,ii) = a(1,1)*b(1,ii);
            c(2,ii) = a(2,1)*b(1,ii);
            c(3,ii) = a(3,1)*b(1,ii);
            c(4,ii) = a(1,1)*b(2,ii);
            c(5,ii) = a(2,1)*b(2,ii);
            c(6,ii) = a(3,1)*b(2,ii);
            c(7,ii) = a(1,1)*b(3,ii);
            c(8,ii) = a(2,1)*b(3,ii);
            c(9,ii) = a(3,1)*b(3,ii);
         end
         return
      end      
   end
   
   
   if na==nb
      if dim == 1
         c = zeros(els,nb);
         for ii=1:nb
            c(1,ii) = a(1,ii)*b(1,ii);
         end
         return
      elseif dim == 2
         c = zeros(els,nb);
         for ii=1:nb
            c(1,ii) = a(1,ii)*b(1,ii);
            c(2,ii) = a(2,ii)*b(1,ii);
            c(3,ii) = a(1,ii)*b(2,ii);
            c(4,ii) = a(2,ii)*b(2,ii);
         end
         return
      else
         c = zeros(els,nb);
         for ii=1:nb
            c(1,ii) = a(1,ii)*b(1,ii);
            c(2,ii) = a(2,ii)*b(1,ii);
            c(3,ii) = a(3,ii)*b(1,ii);
            c(4,ii) = a(1,ii)*b(2,ii);
            c(5,ii) = a(2,ii)*b(2,ii);
            c(6,ii) = a(3,ii)*b(2,ii);
            c(7,ii) = a(1,ii)*b(3,ii);
            c(8,ii) = a(2,ii)*b(3,ii);
            c(9,ii) = a(3,ii)*b(3,ii);
         end
         return
      end      
   end   
   
   c = [];
  

end