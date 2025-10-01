function flag = containsid(this,id_)

    Finger = this.Head;
    while ~isempty(Finger) && ~strcmp(Finger.Data.id,id_)
        Finger = Finger.Next;        
    end % loop until node is found
    flag = ~isempty(Finger);

end % FvmLabEntityList.containsid