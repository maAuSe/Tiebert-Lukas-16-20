function add(this,Data_)

    this.Head = SinglyLinkedListNode(Data_,this.Head);
    this.nElem = this.nElem + 1;
    
end % SinglyLinkedList.add
