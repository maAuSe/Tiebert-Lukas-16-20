function [ir,jc] = calc_ir_jc(nnz,nF,fNbC,n)

   % We reshape fNbC to a two-row vector for our convenience:
   fNbC = reshape(fNbC,2,numel(fNbC)/2);
   % Now fNbC(jF,1) holds nb1; fNbC(jF,2) holds nb2

   % Initialize with monotonic increasing numbers.
   ir = (1:1:nnz)';
   jc = (1:1:nnz)';
	% The first n values, representing the diagonal elements, are now set
   % correctly. We have to adapt the values that are left.
 
    % Loop over faces and write off-diagonal node locations for n1
    % and n2 equations
  cursor = n;
  for jF = 1:nF
      cursor = cursor + 1;                % n1 equation, n2 variable
      ir(cursor) = fNbC(1,jF);
      jc(cursor) = fNbC(2,jF);
            
      cursor = cursor + 1;                % n2 equation, n1 variable
      ir(cursor) = fNbC(2,jF);
      jc(cursor) = fNbC(1,jF);
      
  end
    
     
   ir = cast(ir,'double');
   jc = cast(jc,'double');
    
end