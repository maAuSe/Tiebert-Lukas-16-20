function requirefieldname(var,name,varargin)
        
    if ~isfield(var,name)
        ME = MException('FvmLab:TypeChecker:MissingField',...
           ['\nFvmLab:TypeChecker:MissingField',...
            '\nThe supplied data misses a field with name "%s"'...
            '\n'],name);
        throwAsCaller(ME);
    end
    for iLoop = 1:length(varargin)
        if ~isfield(var,varargin{iLoop})
            ME = MException('FvmLab:TypeChecker:MissingField',...
               ['\nFvmLab:TypeChecker:MissingField',...
                '\nThe supplied data misses a field with name "%s"'...
                '\n'],name);
            throwAsCaller(ME);
        end
    end
        
end