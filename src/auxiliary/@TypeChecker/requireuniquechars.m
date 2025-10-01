function requireuniquechars(var)
% PRE: var is cell array of chars
        
    for iLoop = 1:numel(var)
        testelem = var{iLoop};
        for iLoop2 = iLoop+1:numel(var)
            testelem2 = var{iLoop2};
            if strcmp(testelem,testelem2)
                ME = MException('FvmLab:TypeChecker:DuplicateData',...
                   ['\nFVMLab:TypeChecker:DuplicateData',...
                    '\nThe supplied data elements should not have duplicates'...
                    '\n']);
                throwAsCaller(ME);
            end
        end            
    end
        
end
