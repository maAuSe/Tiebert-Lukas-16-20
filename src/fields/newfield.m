function Obj = newfield(ZoneObj,data)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This factory routine provides a wrapper for easy creation of Field
% objects initialized from the given data. The supplied data is analysed to
% retrieve the rank and dimension. If elcount is supplied, that resulting
% value will have that number of elements.
% POST:
%  A valid Value object iff DATA is of class double and meets one of the
%  following conditions:
%  (1) a = size(data);
%      a(1) == 1 && a(2) == 1
%  (2) a = size(data);
%      (a(1) == (1||2||3) && a(2) == 1) || 
%      (a(2) == (1||2||3) && a(1) == 1)
%  (3) a = size(data)
%      a(1) == a(2) && a(1) == (1||2||3)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   TypeChecker.requireclass(ZoneObj,'IZone')
   TypeChecker.requireclass(data,'double');
   TypeChecker.requirenonempty(data);
   
   dom = getdomain(ZoneObj);
   dim = dom.dim;
   
   ds = size(data);

   
   rank = 0;
   elsize = eldsize(rank,dim);
   if ds(1)==elsize % scalar
      Obj = Field(ZoneObj,rank);
      elcount = getelcount(ZoneObj);
      if ds(2) == 1
         set(Obj,repmat(data,1,elcount));
      else
         set(Obj,data);
      end
      return
   end

   
   rank = 1;
   elsize = eldsize(rank,dim);
   if ds(1)==elsize % vector
      Obj = Field(ZoneObj,rank);
      elcount = getelcount(ZoneObj);
      if ds(2) == 1
         set(Obj,repmat(data,1,elcount));
      else
         set(Obj,data);
      end
      return
   end
   
   
   rank = 2;
   elsize = eldsize(rank,dim);
   if ds(1)==elsize % tensor
      Obj = Field(ZoneObj,rank);
      elcount = getelcount(ZoneObj);
      if ds(2) == 1
         set(Obj,repmat(data,1,elcount));
      else
         set(Obj,data);
      end
      return
   end
   
      
   error('The supplied %gx%g data cannot be fitted into a scalar/vector/tensor value',ds(1),ds(2))
   

end