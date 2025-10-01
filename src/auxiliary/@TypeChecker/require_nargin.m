function require_nargin(val,requiredval)
        
    if val~=requiredval
        error('FVMLAB:InterfaceViolation',...
            'The function was called with the wrong number of input arguments.\n')
    end
        
end