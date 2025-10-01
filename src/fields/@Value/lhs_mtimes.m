function c = lhs_mtimes(a,b)
% lhs_mtimes(a,b) = a * b
% class(a) == 'Value'
% PRE:
%  class(b) == 'Value' || or lower in the hierarchy
% POST
%  class(c) == 'Value'
%

   % Do the operation on the data
   [datac rank dim] = pd_mult(a,b);
   
   % Create the required object
   c = Value(rank,dim);
	set(c,datac);

end