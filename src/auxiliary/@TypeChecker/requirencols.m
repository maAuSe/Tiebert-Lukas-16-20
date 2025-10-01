function requirencols(var,colcount)
        
    if size(var,1) ~= colcount
        ME = MException('FvmLab:TypeChecker:UnexpectedColCount',...
           ['\nFVMLab:TypeChecker:UnexpectedColCount',...
            '\nThe supplied data does not have the required column count of %gx%g'...
            '\n'],colcount);
        throwAsCaller(ME);
    end
        
end