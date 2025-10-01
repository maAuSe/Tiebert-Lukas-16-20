function [d2SLoc,d2S] = inverttopology(nS,nD,s2DLoc,s2D)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% FUNCTION [d2SLoc,d2S] = invertTopology(nS,nD,s2DLoc,s2D)
%
% Calculates the inverse indexing arrays.
%
% Input paremeters are:
%  nS         : number of elements in source
%  nD         : number of elements in destination
%  s2DLoc     : starting location of the destinations
%  s2D        : all destinations
%
% TODO: Validate thoroughly
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   % Expand s2DLoc to d2S
   d2S = cast(zeros(size(s2D)),'int32');
   thisLoc = s2DLoc(1);
   for jS = 1:nS
      nextLoc = s2DLoc(jS+1);
      for jDLoc=thisLoc:nextLoc-1
         d2S(jDLoc) = jS;
      end
      thisLoc = nextLoc;
   end
   
   % Sort s2D and apply sorting to eS2DLoc
   [d2SNdx ndx] = sort(s2D);
   d2S = d2S(ndx);
   
   % Compress d2SNdx to d2SLoc
   d2SLoc = cast(zeros(1,nD+1),'int32');
   lastD = int32(0);
   for iCon = 1:numel(d2SNdx)
      currentD = d2SNdx(iCon);
      if currentD ~= lastD
         d2SLoc(currentD) = iCon;
         lastD = currentD;
      end
   end
   d2SLoc(end) = numel(d2SNdx) + 1;
   
end



   