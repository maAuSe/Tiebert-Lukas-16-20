function setid(This,id)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% setid method for setting the id of This entity
% PRE:
%   1)  "id" must be of char type (empty or filled).
% POST:
%   "This.id" is set to "id" IFF "id" is empty or well formed.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if isempty(id) % Empty id
        This.id = id;
    elseif FvmLabEntity.idiswellformed(id) % Well-formed id
        This.id = id;
    else
        error(['The supplied id ''%s'' must be empty, or must consist of '...
            'a series of valid Matlab variable names seperated by ":"'],id);
    end

end % FvmLabEntity.setid