function Data_ = remove(this,Data_)
% Removes the first occurence of Data_
% Removeall call Remove recursively until nothing is removed

    Finger = this.Head;
    Previous = [];
    while ~isempty(Finger) && ~eq(Finger.Data,Data_)
        Previous = Finger;
        Finger = Finger.Next;
    end % loop until node is found

    if ~isempty(Finger)
    % We found an element to remove
        if isempty(Previous)
            % Remove first element
            this.Head = Finger.Next;
        else
            % Remove other element
            setnext(Previous,Finger.Next)
        end
        this.nElem = this.nElem - 1;
        Data_ = Finger.Data;
    else
        Data_ = [];
    end %

end % SinglyLinkedList.remove