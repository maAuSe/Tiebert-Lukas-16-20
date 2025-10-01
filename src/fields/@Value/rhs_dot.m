function c = rhs_dot(b,a)
% rhs_dot(b,a) = a • b
% class(a) == 'Value'
% PRE:
%  class(b) == 'Value' || or lower in the hierarchy
% POST
%  class(c) == 'Value'
%

   % Do the operation on the data
   [datac rank dim] = pd_dot(a,b);
   
   % Create the required object
   c = Value(rank,dim);
	set(c,datac);

end