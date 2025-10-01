function result = getelementsbyid(this,id_)
% Returns a cellArray containing all elements with given id
    
    result = cell(1,this.nElem);
    nFound = 0;
    Finger = this.Head;
    while ~isempty(Finger)
        if strcmp(Finger.Data.id,id_)
            nFound = nFound + 1;
            result{nFound} = Finger.Data;
        end
        Finger = Finger.Next;        
    end % loop all entries
    if nFound == 0
        result = {};
    else
        result = result(1:nFound);
    end
    

end % FvmLabEntityList.getelementsbyid