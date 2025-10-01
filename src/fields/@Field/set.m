function set(This,data,c1,c2)

   if int32(size(data,2)) ~= This.elcountzone
         error('The supplied data does not match the size required by the zone of this field.')
   end

   if nargin == 2
      c1 = [];
      c2 = [];
   end
   if nargin == 3
      c2 = [];
   end
   
   set@Value(This,data,c1,c2)

end