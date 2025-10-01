classdef SinglyLinkedList < handle
% A linked list with a single fixed iterator the can be used to acces all
% the elements.

properties(SetAccess=protected,GetAccess=public)
    nElem    
end % PROPERTIES


properties(SetAccess=protected,GetAccess=protected)
    % We want to avoid access to SinglyLinkedListNode objects by the user.
    % The user only manipulated the data!
    Head
    
end % PROTECTED PROPERTIES



methods
    
    function this = SinglyLinkedList()
    
        this.Head = [];
        this.nElem = 0;
        
	end % SinglyLinkedList CONSTRUCTOR
    
    
    flag    = isempty(this)    
    clearlist(this)
    add(this,Data_)
    Data_   = removefromhead(this)
    Data_   = peek(this)
    flag    = contains(this,Data_)
    Data_   = remove(this,Data_)
    removeall(this,Data_)
    it      = elements(this)
       
end % METHODS



end % classdef SinglyLinkedList