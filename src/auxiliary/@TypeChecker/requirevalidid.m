function requirevalidid(id_)
        
    if ~FvmLabEntity.idiswellformed(id_)
        error('FVMLab:TypeChecker:WrongIDType',...
            'The supplied id_ is not a valid FvmLabID');
    end
        
end