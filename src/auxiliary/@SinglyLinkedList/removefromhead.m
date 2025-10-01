function Data_ = removefromhead(this)

    temp = this.Head;
    if ~isempty(temp)
        this.Head = temp.Next;
        this.nElem = this.nElem - 1;
        Data_ = temp.Data;
    else
        Data_ = [];
    end
    
end % SinglyLinkedList.removefromhead