function requireclass(var,classString)
        
    if ~isa(var,classString)
        ME = MException('FVMLab:TypeChecker:WrongClass',...
            ['\nFVMLab:TypeChecker:WrongClass',...
             '\nThe supplied class "%s" is different from the expected class "%s"'...
             '\n'],class(var),classString);
        throwAsCaller(ME);
    end
        
end