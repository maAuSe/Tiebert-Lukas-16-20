function setdescr(This,descr_)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% setdescr method for setting the description of This entity
% PRE:
%   /
% POST:
%   "This.descr_" is set to "descr_" IFF "descr_" is a character array.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if isa(descr_,'char') % Description is charater array
        This.descr = descr_;
    else
        error('The supplied argument must be a character array');
    end

    This.descr = ''; % Initialize with no description

end % FvmLabEntity.setdescr