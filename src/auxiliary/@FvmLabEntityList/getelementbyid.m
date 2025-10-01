function result = getelementbyid(this,id_)

    Finger = this.Head;
    while ~isempty(Finger) && ~strcmp(Finger.Data.id,id_)
        Finger = Finger.Next;        
    end % loop until node is found
    if isempty(Finger)
        result = [];
    else
        result = Finger.Data;
    end

end % FvmLabEntityList.getelementbyid