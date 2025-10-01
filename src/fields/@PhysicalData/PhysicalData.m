classdef PhysicalData < handle



properties (GetAccess=public,SetAccess=protected)

   rank        % int32
   % Rank of data (0 = scalar, 1=vector, 2=tensor)
   dim         % int32
   % Physical dimension of data (1 ~> 1Dv, 2 ~> (2Dv or 2Dt) ~> 3=(3Dv or 3Dt) 
   
   elsize      % int32
   % = eldsize(rank,dim).
   % Holds the number of double precision memory locations to store a single element
   
   isscal      % logical
   % = (rank == 0)
   isvec       % logical
   % = (rank == 1)
   istens      % logical
   % = (rank == 2)
    
end % PROPERTIES




methods (Access=public)

    
   function Obj = PhysicalData(rank,dim)
   %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   % PhysicalData constructor.
   % PRE:
   %   rank and data can be anything
   % POST:
   %   A valid PhysicalData object
   %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      if rank == 0
         Obj.rank = int32(0);
         Obj.dim = setdim(dim);
         Obj.elsize = int32(dim^0);
         Obj.isscal = true;
         Obj.isvec = false;
         Obj.istens = false;         
      elseif rank == 1
         Obj.rank = int32(1);
         Obj.dim = setdim(dim);
         Obj.elsize = int32(dim^1);
         Obj.isscal = false;
         Obj.isvec = true;
         Obj.istens = false;         
      elseif rank == 2
         Obj.rank = int32(2);
         Obj.dim = setdim(dim);
         Obj.elsize = int32(dim^2);
         Obj.isscal = false;
         Obj.isvec = false;
         Obj.istens = true;
      else
         error('Rank must be 0,1 or 2')
      end
      
            
      function sdim = setdim(dim)
      % Subfunction to validate dimension and return correct class
         if dim == 1
            sdim = int32(1);
         elseif dim == 2
            sdim = int32(2);
         elseif dim == 3
            sdim = int32(3);
         else
            error('The supplied dimension must be 1, 2 or 3')
         end
      end

   end % PhysicalData constructor
       
   flag = isscalar(a);
   flag = isvector(a);
   flag = istensor(a);

end % PUBLIC METHODS




methods (Access=public,Abstract=true)

   elcount = getelcount(This)
   haselflag = haselements(This)

end % ABSTRACT METHODS




methods (Static=true)
    
     
    function err_binop_unexpected_lhs(op,lhsClass,rhsClass)
        
        ME = MException('FvmLab:PhysicalData:BinOp:UnexpectedInput',...
               ['FvmLab:PhysicalData:BinOp:UnexpectedInput',...
                '\n\nThe operation ''%s'' is not supported for the supplied left hand side A.',...
                '\nThe classes of the supplied data were:\n',...
                ' A --> ''%s''\n',...
                ' B --> ''%s''\n',...
                '\nRefer to the manual for valid use of this operation.',...
                '\n'
               ],op,lhsClass,rhsClass);
        throwAsCaller(ME);
        
    end % err_binop_unexpected_lhs
        
    
    function err_binop_unexpected_rhs(op,lhsClass,rhsClass)
        
        ME = MException('FvmLab:PhysicalData:BinOp:UnexpectedInput',...
               ['FvmLab:PhysicalData:BinOp:UnexpectedInput',...
                '\n\nThe operation ''%s'' is not supported for the supplied right hand side B.',...
                '\nThe classes of the supplied data were:\n',...
                ' A --> ''%s''\n',...
                ' B --> ''%s''\n',...
                '\nRefer to the manual for valid use of this operation.',...
                '\n'
               ],op,lhsClass,rhsClass);
        throwAsCaller(ME);
        
    end % err_binop_unexpected_rhs
    

    function err_binop_unequal_rank(op,lhsRank,rhsRank)
        
        ME = MException('FvmLab:PhysicalData:BinOp:UnequalRank',...
               ['FvmLab:PhysicalData:BinOp:UnequalRank',...
                '\n\nThe operation ''%s'' is only supported for equal rank operands.',...
                '\nThe ranks of the supplied data were:\n',...
                ' A --> %g\n',...
                ' B --> %g\n',...
                '\nRefer to the manual for valid use of this operation.',...
                '\n'
               ],op,lhsRank,rhsRank);
        throwAsCaller(ME);
        
    end % err_binop_unequal_rank
            
    
    function err_binop_unequal_dim(op,lhsDim,rhsDim)
        
        ME = MException('FvmLab:PhysicalData:BinOp:UnequalDim',...
               ['FvmLab:PhysicalData:BinOp:UnequalDim',...
                '\n\nThe operation ''%s'' is only supported for equal spacial dimension operands.',...
                '\nThe spacial dimensions of the supplied data were:\n',...
                ' A --> %g\n',...
                ' B --> %g\n',...
                '\nRefer to the manual for valid use of this operation.',...
                '\n'
               ],op,lhsDim,rhsDim);
        throwAsCaller(ME);
        
    end % err_binop_unequal_dim
    
    
    function err_binop_unequal_unit(op,lhsUnitString,rhsUnitString)
        
        ME = MException('FvmLab:PhysicalData:BinOp:UnequalUnit',...
               ['FvmLab:PhysicalData:BinOp:UnequalUnit',...
                '\n\nThe operation ''%s'' is only supported for equal units of operands.',...
                '\nThe units of the supplied data were:\n',...
                ' A --> %g\n',...
                ' B --> %g\n',...
                '\nRefer to the manual for valid use of this operation.',...
                '\n'
               ],op,lhsUnitString,rhsUnitString);
        throwAsCaller(ME);
        
    end % err_binop_unequal_unit
    
    
    function err_binop_unequal_zone(op,lhsZoneString,rhsZoneString)
        
        ME = MException('FvmLab:PhysicalData:BinOp:UnequalZone',...
               ['FvmLab:PhysicalData:BinOp:UnequalZone',...
                '\n\nThe operation ''%s'' is only supported for equal zone operands.',...
                '\nThe zones of the supplied data were:\n',...
                ' A --> ''%s''\n',...
                ' B --> ''%s''\n',...
                '\nRefer to the manual for valid use of this operation.',...
                '\n'
               ],op,lhsZoneString,rhsZoneString);
        throwAsCaller(ME);
        
    end % err_binop_unequal_zone
    
    
    function err_set_unexpected_input(inputClass,expectedClass)
        
        ME = MException('FvmLab:PhysicalData:Set:UnexpectedInput',...
               ['FvmLab:PhysicalData:Set:UnexpectedInput',...
                '\n\nThe supplied data must be of type ''%s''.',...
                '\nThe class of the supplied data was ''%s''',...
                '\n\nRefer to the manual for valid use of this method.',...
                '\n'
               ],expectedClass,inputClass);
        throwAsCaller(ME);
        
    end % err_set_unexpected_input
    
    
    function err_set_incorrect_datasize(suppliedNCols,requiredNCols)
        
        ME = MException('FvmLab:PhysicalData:Set:IncorrectDataSize',...
               ['FvmLab:PhysicalData:Set:IncorrectDataSize',...
                '\n\nThe number of columns in the supplied data must match the requirements of the object.',...
                '\nThe required number of columns is %d',...
                '\nThe supplied number of columns is %d',...                
                '\n\nRefer to the manual for valid use of this method.',...
                '\n'
               ],requiredNCols,suppliedNCols);
        throwAsCaller(ME);
        
    end % err_set_unexpected_input
    
    
    function err_set_incorrect_dim(suppliedDim,requiredDim)
        
        ME = MException('FvmLab:PhysicalData:Set:IncorrectDim',...
               ['FvmLab:PhysicalData:Set:IncorrectDim',...
                '\n\nThe number of rows in the supplied data must match the requirements of the object.',...
                '\nThe required number of rows is %d',...
                '\nThe supplied number of rows is %d',...
                '\n\nRefer to the manual for valid use of this method.',...
                '\n'
               ],requiredDim,suppliedDim);
        throwAsCaller(ME);
        
    end % err_set_unexpected_input
    
    function err_unsupported_op(reasonStr)
        
        ME = MException('FvmLab:PhysicalData:UnsupportedOperation',reasonStr);
        throwAsCaller(ME);
        
    end
    
    % Unsupported operation error


end % STATIC METHODS



end












