function c = rhs_plus(b,a)
% rhs_plus(b,a) = a + b
%  class(b) == 'Field'
% PRE:
%  class(a) == 'Field' || or lower in the hierarchy
% POST
%  class(c) == 'Field'
%

   zone = getzone(b);
   
   if isa(a,'IField') && zone ~= getzone(a)
      error('Attempting to add fields on different zones.')
   end

   % Do the operation on the data
   [datac rank dim] = pd_add(a,b);
   
   % Create the required object
   c = Field(zone,rank);
	set(c,datac);

end