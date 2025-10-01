function Data_ = removeid(this,id_)
% Removes the first FvmLabEntity with .id equal to id_

    Finger = this.Head;
    Previous = [];
    while ~isempty(Finger) && ~strcmp(Finger.Data.id,id_)
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

end % FvmLabEntityList.removeid