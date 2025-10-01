function Obj = fromfile(pathToFile)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% STATIC METHOD FvMesh.fromfile(pathToFile)
% 
% Facory method for creating an FvMesh object from data in a .mat file.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   fileContent = load(pathToFile,'meshStruct');
   Obj = FVMLabMeshObj.fromstruct(fileContent.mStruct,'fill');
   clear('fileContent');

end