function Data_ = nextelement(this)
% No error checks are made (by purpose => speed)
% Use hasmoreelements to check this beforehand!

    Data_ = this.Current.Data;    
    this.Current = this.Current.Next;    
    
end % SinglyLinkedListIterator.nextelement