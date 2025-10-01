classdef SinglyLinkedListNode < handle
 
properties(SetAccess=public,GetAccess=public)
    Data
end % PUBLIC PROPERTIES
    
    
    
properties(SetAccess=protected,GetAccess=public)
    Next
end % PROPERTIES    


methods

    function this = SinglyLinkedListNode(Data_,NextNode_)
        
    	this.Data = Data_;
        if nargin == 2
            % not empty: node is element followed by next.
            if ~isempty(NextNode_) && ~isa(NextNode_,'SinglyLinkedListNode')
                error('Supplied data must be a SinglyLinkedListNode');
            end
            this.Next = NextNode_; 
        else
            % empty: node is tail node
            this.Next = []; 
        end
        
    end % SinglyLinkedListNode CONSTRUCTOR
    
    % NOTE: the setnext method was forgotten inthe spec
    setnext(this,NextNode_)

end % METHODS
    
    
end % classdef SinglyLinkedListNode