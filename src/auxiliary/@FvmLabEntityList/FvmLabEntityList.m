classdef FvmLabEntityList < SinglyLinkedList
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% CLASS FvmLabEntityList
%
% This class adds to a SinglyLinkedList a convenient mechanism for storing
% FvmLabEntity objects and retrieving them by their Id. Items with empty Id
% are not allowed.
% When the "noDuplicatesFlag" is on, no duplicates are allowed in the
% strong sense that adding items with duplicated Id's is illegal.
% 
%
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


properties
    
    noDuplicatesFlag    % logical(1,1)
    % When this flag is on, the list is searched to see if the element to
    % be added is already in the list.
    
end % PROPERTIES



methods (Access=public)
    
    function this = FvmLabEntityList()
    
        this = this@SinglyLinkedList();
        setnoduplicatesflag(this,true);
        
	end % FvmLabEntityList CONSTRUCTOR
        
    % Overloaded methods
    add(this,Data_)
    
    % Added functionality
    setnodusplicatesflag(this,flag)
    result = getelementbyid(this,id_)
    result = getelementsbyid(this,id_)
    flag   = containsid(this,id_)
    Data_  = removeid(this,id_)
       
end % METHODS



end % classdef FvmLabEntityList