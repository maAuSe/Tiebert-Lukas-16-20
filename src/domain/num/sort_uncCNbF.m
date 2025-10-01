function cNbF = sort_uncCNbF(nPc,cNbFLoc,cNbF,fNbVLoc,fNbV)
% This function sorts the cNbF array in such a way that consecutive faces
% in the list of each cell are connected in the mesh. The faces might get
% sorted clockwise or counter-clockwise for different cells.

    
   fNbV = reshape(fNbV,fNbVLoc,[]);

   for iCell = 1:1:nPc
      startPtr = cNbFLoc(iCell);
      endPtr = cNbFLoc(iCell + 1) - 1;
      lastSortedPtr = startPtr;
      lastSortedFace = cNbF(lastSortedPtr);
      currCheckFacePtr = lastSortedPtr;
      lastSortedV1 = fNbV(1,lastSortedFace);
      lastSortedV2 = fNbV(2,lastSortedFace);
      while lastSortedPtr < endPtr
         % move current to next position
         currCheckFacePtr = currCheckFacePtr + 1;
         currCheckFace = cNbF(currCheckFacePtr);
         currCheckV1 = fNbV(1,currCheckFace);
         currCheckV2 = fNbV(2,currCheckFace);
         % check if current is connected with lastSorted
         if lastSortedV1==currCheckV1 || lastSortedV1==currCheckV2 || ...
               lastSortedV2==currCheckV1 || lastSortedV2==currCheckV2
            % connected face found:
            lastSortedV1 = currCheckV1;
            lastSortedV2 = currCheckV2;
            lastSortedFace = currCheckFace;
            lastSortedPtr = lastSortedPtr + 1;
            if lastSortedPtr ~= currCheckFacePtr
               % The list is not sorted, so entries must be swapped
               cNbF(currCheckFacePtr) = cNbF(lastSortedPtr);
               currCheckFacePtr = lastSortedPtr;
               % (swap thisFace with the face residing at lastSortedPtr location)
            end %
            cNbF(lastSortedPtr) = lastSortedFace;
         else % no connection found
            % nothing should happen?
         end
         % if not so: thisFace
      end % loop over connected faces
   end % loop over cells

end