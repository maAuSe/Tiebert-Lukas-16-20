function Data_ = peek(this)
    
    temp = this.Head;
    if ~isempty(temp)
        Data_  = this.Head.Data;        
    else
        Data_ = [];
    end
    
end % SinglyLinkedList.peek