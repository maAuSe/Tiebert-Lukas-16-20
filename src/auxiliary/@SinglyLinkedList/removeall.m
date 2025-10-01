function removeall(this,Data_)
% PRE: algorithm does not work with isempty(Data_)
%   TODO: in that case switch to private removeallempty method. Reason for
%   this hack is that matlab forces us to use empty arrays in place of NULL 
%   references.

    if isempty(Data_)
        error('This method cannot remove all empty data.')
    end
    while ~isempty(Data_)
        Data_ = remove(this,Data_);
    end

end % SinglyLinkedList.removeall