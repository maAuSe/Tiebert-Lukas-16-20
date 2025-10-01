function c = lhs_plus(a,b)
% lhs_plus(a,b) = a + b
%  class(a) == 'Field'
% PRE:
%  class(b) == 'Field' || or lower in the hierarchy
% POST
%  class(c) == 'Field'
%

   zone = getzone(a);

   if isa(b,'IField') && zone ~= getzone(b)
      error('Attempting to add fields on different zones.')
   end

   % Do the operation on the data
   [datac rank dim] = pd_add(a,b);
   
   % Create the required object
   c = Field(zone,rank);
	set(c,datac);

end