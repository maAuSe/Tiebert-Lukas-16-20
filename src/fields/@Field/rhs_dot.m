function c = rhs_dot(b,a)
% rhs_dot(b,a) = a • b
% class(a) == 'Field'
% PRE:
%  class(b) == 'Field' || or lower in the hierarchy
% POST
%  class(c) == 'Field'
%

   zone = getzone(b);
   
   if isa(a,'IField') && zone ~= getzone(a)
      error('Attempting to dot fields on different zones.')
   end

   % Do the operation on the data
   [datac rank dim] = pd_dot(a,b);
   
   % Create the required object
   c = Field(zone,rank);
	set(c,datac);

end