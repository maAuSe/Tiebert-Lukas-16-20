function flag = idiswellformed(id_)
% TODO: rename to validate id without the flag
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% PRE:
%   1)  "id_" must be of char type
% POST:
%   1)  "flag" indicates wether the "id_" is well-formed in the context
%       of FvmLab.
% NOTE: Well-formedness rules for ID's are the same as for Matlab
% variables, however with the exception that colons ':' are allowed as
% seperator to allow for nested id's that belong to some group.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    flag = true;
    while true
        [str, id_] = strtok(id_, ':'); %#ok<STTOK>
        if isempty(str),  break;
        elseif ~isvarname(str)
            flag = false;
            break;
        end            
    end

end % FvmLabEntity.idiswellformed