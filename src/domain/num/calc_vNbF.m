function [vNbFLoc vNbF] = calc_vNbF(nV,nF,fNbVLoc,fNbV)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% FUNCTION [vNbFLoc vNbF] = calc_vNbF(fNbVLoc,fNbV,nF,nV)
%
% Calculates the topology information for indexing from vertices to faces
% on an unstructured mesh. Input paremeters are:
%  nV          : number of vertices
%  nF          : number of faces
%  fNbVLoc     : starting locations of neighbors in fNbV
%  fNbV        : contains the vertex neighbors of all faces
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    vNbF = repmat(int32(1:1:nF),2,1);
    vNbF = reshape(vNbF,1,[]);

    [vNbFNdx ndx] = sort(fNbV);
    vNbF = vNbF(ndx);

    % Compress (vNbFNdx    ->    vNbF) 
    %   into   (vNbF       ->    vNbF)
    vNbFLoc = zeros(1,nV+1);
    lastFace = int32(0);
    for k = int32(1):length(vNbF)
        currentFace = vNbFNdx(k);
        if currentFace ~= lastFace
            vNbFLoc(currentFace) = k;
            lastFace = currentFace;
        end
    end
    vNbFLoc(end) = k+1; % Last element

end



