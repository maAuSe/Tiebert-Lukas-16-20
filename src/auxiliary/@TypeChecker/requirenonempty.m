function requirenonempty(var)
        
    if isempty(var)
        ME = MException('FvmLab:TypeChecker:UnexpectedEmptyData',...
           ['\nFVMLab:TypeChecker:UnexpectedEmptyData',...
            '\nThe supplied data cannot be empty'...
            '\n']);
        throwAsCaller(ME);
    end
        
end