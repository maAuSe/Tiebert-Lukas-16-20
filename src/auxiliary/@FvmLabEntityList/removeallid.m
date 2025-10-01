function removeallid(this,id_)

    Data_ = removeid(this,id_);
    while ~isempty(Data_)
        id_ = Data_.id_;
        Data_ = removeid(this,id_);
    end

end % FvmLabEntityList.removeallid