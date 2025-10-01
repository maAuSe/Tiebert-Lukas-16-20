function tofile(This,pathToFile)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% PUBLIC METHOD tofile(FvMeshObj,pathToFile)
% 
% Saves an existing FvMesh object to a .mat file with full path string
% specification. 
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   mStruct = struct();

   mStruct.dim    = This.dim;
   mStruct.nV     = This.nV;
   mStruct.nF     = This.nF;
   mStruct.nC     = This.nC;
   mStruct.nIv    = This.nI;
   mStruct.nBv    = This.nBv;
   mStruct.nIf    = This.nIf;
   mStruct.nBf    = This.nBf;
   
   mStruct.vCoord       = This.vCoord;

   mStruct.fNbVLoc      = This.fNbVLoc;
   mStruct.fNbV         = This.fNbV;
   mStruct.fNbCLoc      = This.fNbCLoc;
   mStruct.fNbC         = This.fNbC;
   mStruct.maxNoCNbF    = This.maxNoCNbF;
   mStruct.maxNoCNbV    = This.maxNoCNbV;   
   
   mStruct.noNamedBfaceZones ...
                        = This.noNamedBfaceZones;
	mStruct.bfaceData    = This.bfaceData;

   % Leave out '-v6' to use compression
   save(pathToFile,'-v6','mStruct');

end