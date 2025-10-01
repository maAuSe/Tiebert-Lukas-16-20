classdef FvmLabEntity < handle
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% CLASS FvDomain 
%
% This class defines entities with an id and with a description.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


properties (GetAccess=public,SetAccess=protected)
    
    id                  % char
    %   Well-formed identifier for this entity.    
    descr               % char
    %   Description of this entity.

end % PROTECTED PROPERTIES




methods (Access=public)
    
    function This = FvmLabEntity(id)
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % FvmLabEntity constructor
    % PRE:
    %   1)  "id" must be of char type (empty or filled).
    % POST:
    %   FvmLabEntity "obj" is constructed IFF "id_" is empty or well
    %   formed. 
	%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        setid(This,id);
        setdescr(This,char([]));        
        
    end % FvmLabEntity CONSTRUCTOR
    
    setid(This,id)
    id = getid(This)
    setdescr(This,descr)
        
end % PUBLIC METHODS




methods (Static)
    
    flag = idiswellformed(id_)

end % STATIC METHODS



end % FvmLabEntity CLASSDEF