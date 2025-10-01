classdef TypeChecker
    

    methods (Static)
    
        requireclass(var,classString)
        
        requiresingleelem(var)
        
        requirerowarray(var,rowlength)
        
        requirecolarray(var,collength)
        
        requirearraysize(var,rowlength,collength)
        
        requirenrows(var,rowcount)
        
        requirencols(var,colcount)
        
        requirevalidid(id_)
        
        requirefieldname(var,name,varargin)
                
        requirenonempty(var)        
        
        requireuniquechars(var)
        
        require_nargin(val,requiredval)
        
        require_inputclass(var,requiredclass)
        
        require_rank(var,allowedranks)
        
        require_validrange(range,varargin)
    
    end % STATIC
    
    
    
end % classdef TypeChecker