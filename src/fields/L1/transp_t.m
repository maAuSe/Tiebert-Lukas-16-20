function a = transp_t(a,n,dim)
% A is a tensor array of dimension DIM stored in column major.
% N is the number of tensors in the array

if dim == 1
   return   
end


if dim == 2
   tmp = 0;
   for ii=1:n
      tmp = a(2);
      a(2) = a(3);
      a(3) = tmp;
   end
   return
end


if dim == 3
   tmp = 0;
   for ii=1:n
      tmp = a(2);
      a(2) = a(4);
      a(4) = tmp;
      tmp = a(3);
      a(3) = a(7);
      a(7) = tmp;
      tmp = a(6);
      a(6) = a(8);
      a(8) = tmp;
   end
   return
end




end