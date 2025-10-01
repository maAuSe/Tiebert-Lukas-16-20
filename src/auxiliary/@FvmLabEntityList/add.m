function add(this,Data_)

    if ~isa(Data_,'FvmLabEntity')
        error('A FvmLabList can only contain ''FvmLabEntity'' objects.');
    end
    if this.noDuplicatesFlag
        if containsid(this,Data_.id)
            % || contains(this,Data_) : This check is unnecessary:
            % identical objects always have identical id's!
            error('No duplicates are allowed when the noDuplicatesFlag is on.')
        end
    end
    if isempty(Data_.id)
        error('No empty entries are allowed in a FvmLabEntityList.')
    end
    add@SinglyLinkedList(this,Data_)

end % FvmLabEntityList.add
