function requirenrows(var,rowcount)
        
    if size(var,1) ~= rowcount
        ME = MException('FvmLab:TypeChecker:UnexpectedRowCount',...
           ['\nFVMLab:TypeChecker:UnexpectedRowCount',...
            '\nThe supplied data does not have the required row count of %g'...
            '\n'],rowcount);
        throwAsCaller(ME);
    end
        
end