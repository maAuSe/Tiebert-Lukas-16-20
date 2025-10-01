function m = maxabs(a)

   elsize = size(a,1);
   elcount = size(a,2);
   
   m = 0;
   	
   for jc = 1:elcount
      for ir = 1:elsize
         m = max(m,abs(a(ir,jc)));
      end      
   end

end