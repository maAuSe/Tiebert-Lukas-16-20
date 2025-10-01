function c = lhs_minus(a,b)
% lhs_minus(a,b) = a - b
% class(a) == 'Value'
% PRE:
%  class(b) == 'Value' || or lower in the hierarchy
% POST
%  class(c) == 'Value'
%

   % Do the operation on the data
   [datac rank dim] = pd_subtr(a,b);
   
   % Create the required object
   c = Value(rank,dim);
	set(c,datac);

end