function requiresingleelem(var)
        
    if ~isscalar(var)
        error('FVMLab:TypeChecker:NoSingleValue',...
            'The supplied data is not scalar');
    end
        
end