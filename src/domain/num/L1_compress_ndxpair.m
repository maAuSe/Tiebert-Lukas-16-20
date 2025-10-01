function [ndxLoc comprNdxVec] = L1_compress_ndxpair(ndxLoc,uncomprNdxVec,comprNdxVec)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% FUNCTION [ndxLoc newNdxVec] = L1_compress_ndxpair(
%                 ndxLoc,oldNdxVec,newNdxVec)
%
% Compresses an uncompressed indexing vector (typically constructed from an
% estimate of the maximum connectivity) into a compressed one, that is
% sized according to the actual connectivity. The location vector is
% altered in place. Zero entries of the uncompressed vector are cancelled.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ptrU = ndxLoc(1);    % U suffix stands for *U*ncompressed
ptrC = ptrU;         % C suffix stands for *C*ompressed
for jLoc = 1:numel(ndxLoc)-1
   ptrStartNext = ndxLoc(jLoc+1);
   
   % Loop all entries:
   % if not zero, copy into comprNdxVec and move this ptrStartC
   % if zero, discard
   for jj = ptrU:1:ptrStartNext-1
      tmp = uncomprNdxVec(jj);
      if tmp~=0
         comprNdxVec(ptrC) = uncomprNdxVec(jj);
         ptrC = ptrC+1;
      end      
   end
   
   % Change the target pointer to the location on ndxLoc
   ndxLoc(jLoc+1) = ptrC;
   
   % Move pointer
   ptrU = ptrStartNext;
end