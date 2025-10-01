function requirearraysize(var,rowlength,collength)
        
    if size(var,1) ~= rowlength || size(var,2) ~= collength        
        ME = MException('FvmLab:TypeChecker:UnexpectedArraySize',...
           ['\nFVMLab:TypeChecker:UnexpectedArraySize',...
            '\nThe supplied data does not have the required array size of %gx%g'...
            '\n'],rowlength,collength);
        throwAsCaller(ME);
    end
        
end