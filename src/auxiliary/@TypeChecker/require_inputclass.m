function require_inputclass(var,requiredclass)
        
   if ~strcmp(class(var),requiredclass)
      error('FVMLAB:InterfaceViolation',...
            'The input argument "%s" has a wrong class.\n',inputname(1));
   end
        
end