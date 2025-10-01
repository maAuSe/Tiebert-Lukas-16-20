function c = rhs_minus(b,a)
% rhs_minus(b,a) = a - b
% class(b) == 'Value'
% PRE:
%  class(a) == 'Value' || or lower in the hierarchy
% POST
%  class(c) == 'Value'
%

   % Do the operation on the data
   [datac rank dim] = pd_subtr(a,b);
   
   % Create the required object
   c = Value(rank,dim);
	set(c,datac);

end