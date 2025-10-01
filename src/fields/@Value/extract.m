function data = extract(This)

   % if ~haselements(This)
   %    data = [];
   %    return
   % end

   data = This.data;
	This.data = [];
   This.elcount = 0;
   return

end