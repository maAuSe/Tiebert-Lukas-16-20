classdef FvMesh
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% CLASS FvMesh
%
% A FvMesh holds the minimal set of geometry and mesh data that can be
% processed for creating a domain. A mesh can be created from a struct,
% loaded from file, and saved to file.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

properties (GetAccess=public,SetAccess=public)

   % VALUE PROPERTIES -----------------------------------------------------
   dim      % int32(1,1)
   % Spatial dimension of mesh
   nV       % int32(1,1)
   % Number of Vertices  
   nF       % int32(1,1)
   % Number of Faces
   nPc      % int32(1,1)
   % Number of Cells
   
   nIv      % int32(1,1)
   % Number of Internal Vertices
   nBv      % int32(1,1)
   % Number of Boundary Vertices
   nIf      % int32(1,1)
   % Number of Internal Faces
   nBf      % int32(1,1)
   % Number of Boundary Faces

   
   % GEOMETRY -------------------------------------------------------------
   vCoord      	% double(2,noVertices)
   % Vertices X-coordinates in 1st row
   % Vertices Y-coordinates in 2nd row

   
   % TOPOLOGY -------------------------------------------------------------
   % Note that a single value for the Loc arrays indicates that there is
   % a fixed number of neighbors of the specified kind. The value itself
   % indicates This fixed number.
   fNbVLoc        % int32(1,1)
   % Holds the value 2 for a 2D mesh since every face connects to two
   % vertices.
   fNbV           % int32(1,fNbVLoc*noFaces)
   % Holds numbers of the vertices connected to each face.
   fNbCLoc        % int32(1,1)
   % Always holds the value 2 since a can connect at maximum to two cells.
   fNbC           % int32(1,fNbCLoc*nFaces)
   % Internal faces always connect to two cells. Boundary faces connect to
   % a single cell. The second connection is set to 0. 
   maxNoCNbF      % int32(1,1)
   % Maximum number of cell neighbor faces. Used for allocation of
   % intermediate data structures.
   maxNoCNbV      % int32(1,1)
   % Maximum number of cell neighbor vertices. Used for allocation of
   % intermediate data structures.
   
   % ORDERING -------------------------------------------------------------
   % The order of faces and cells is such that internal entities come first
   % and boundary entities come last. This is an additional specification
   % (on top of the topology) that enables efficient looping over internal 
   % regions and (parts of) boundary regions. Formally it is as follows:
   %   If a face has no cell neighbor
   %   Then
   %     (1) The face number must be at least noIntrnlFaces+1
   %     (2) The vertex numbers of the defining vertices must be at least
   %         noIntrnlVertices+1
   % Note that at This time, the specification has no special ordering for
   % cells lying against the boundary. This implies that looping over these
   % cells cannot be done directly. BUT, it can be done relatively
   % efficient by looping over the boundary faces and retrieving the owner
   % cell.
    

  	% NAMED BOUNDARY FACE ZONES -------------------------------------------    
   noNamedBfaceZones    % int32(1,1)
   % Number of boundary face zones that are defined in This mesh
   bfaceData            % struct(1,noNamedBfaceZones)
   %   .id     : charArray
   %       A valid id specification
   %   .range  : int32(1,2)
   %       Specifies the face range that are named with the given id


end % PROPERTIES




methods (Access=public)

    display(This)    
    tofile(This,pathToFile)
    
end % METHODS




methods(Access=public,Static=true)

    Obj = fromstruct(mStruct)
    Obj = fromfile(pathToFile)
    
end % STATIC METHODS




methods (Access=protected)
    
	function This = FvMesh(mStruct)
	% Protected constructor. For internal use only. Use the static factory
	% methods for creating FvMesh objects.
   
      try         
         TypeChecker.requireclass(mStruct.dim,'int32');
         TypeChecker.requiresingleelem(mStruct.dim);         
         This.dim          = mStruct.dim;

         TypeChecker.requireclass(mStruct.nV,'int32');
         TypeChecker.requiresingleelem(mStruct.nV);
         This.nV  = mStruct.nV;

         TypeChecker.requireclass(mStruct.nF,'int32');
         TypeChecker.requiresingleelem(mStruct.nF);
         This.nF  = mStruct.nF;

         TypeChecker.requireclass(mStruct.nPc,'int32');
         TypeChecker.requiresingleelem(mStruct.nPc);
         This.nPc = mStruct.nPc;

         TypeChecker.requireclass(mStruct.nIv,'int32');
         TypeChecker.requiresingleelem(mStruct.nIv);
         This.nIv = mStruct.nIv;

         TypeChecker.requireclass(mStruct.nBv,'int32');
         TypeChecker.requiresingleelem(mStruct.nBv);
         This.nBv = mStruct.nBv;

         TypeChecker.requireclass(mStruct.nIf,'int32');
         TypeChecker.requiresingleelem(mStruct.nIf);
         This.nIf = mStruct.nIf;

         TypeChecker.requireclass(mStruct.nBf,'int32');
         TypeChecker.requiresingleelem(mStruct.nBf);
         This.nBf = mStruct.nBf;

         
         TypeChecker.requireclass(mStruct.vCoord,'double');
         TypeChecker.requirearraysize(mStruct.vCoord,...
            mStruct.dim,mStruct.nV);
         This.vCoord       = mStruct.vCoord;
         
         
         TypeChecker.requireclass(mStruct.fNbVLoc,'int32');
         TypeChecker.requiresingleelem(mStruct.fNbVLoc);
         This.fNbVLoc      = mStruct.fNbVLoc;
         
         TypeChecker.requireclass(mStruct.fNbV,'int32');
         TypeChecker.requirerowarray(mStruct.fNbV,...
            mStruct.fNbVLoc*mStruct.nF);
         This.fNbV         = mStruct.fNbV;
   
         TypeChecker.requireclass(mStruct.fNbCLoc,'int32');
         TypeChecker.requiresingleelem(mStruct.fNbCLoc);
         This.fNbCLoc      = mStruct.fNbCLoc;
         
         TypeChecker.requireclass(mStruct.fNbC,'int32');
         TypeChecker.requirerowarray(mStruct.fNbC,...
            mStruct.fNbCLoc*mStruct.nF);
         This.fNbC         = mStruct.fNbC;
         
         TypeChecker.requireclass(mStruct.maxNoCNbF,'int32');
         TypeChecker.requiresingleelem(mStruct.maxNoCNbF);
         This.maxNoCNbF    = mStruct.maxNoCNbF;

         TypeChecker.requireclass(mStruct.maxNoCNbV,'int32');
         TypeChecker.requiresingleelem(mStruct.maxNoCNbV);
         This.maxNoCNbV    = mStruct.maxNoCNbV;

         TypeChecker.requireclass(mStruct.noNamedBfaceZones,'int32');
         TypeChecker.requiresingleelem(mStruct.noNamedBfaceZones);
         This.noNamedBfaceZones  = mStruct.noNamedBfaceZones;
         
         if mStruct.noNamedBfaceZones == 0
            % Create empty struct
            This.bfaceData = struct('id',{},'range',{});
         else
            % Validate existing struct of bfaceData
            TypeChecker.requirefieldname(mStruct.bfaceData,...
               'id','range');
            TypeChecker.requirearraysize(mStruct.bfaceData,...
               1,mStruct.noNamedBfaceZones);
            for iLoop = 1:mStruct.noNamedBfaceZones
               TypeChecker.requirevalidid(...
                  mStruct.bfaceData(iLoop).id);
               
               TypeChecker.requireclass(...
                  mStruct.bfaceData(iLoop).range,'int32');
               TypeChecker.requirenonempty(...
                  mStruct.bfaceData(iLoop).range);
               TypeChecker.requirearraysize(...
                  mStruct.bfaceData(iLoop).range,2,1);
               TypeChecker.require_validrange(...
                  mStruct.bfaceData(iLoop).range,...
                  mStruct.nIf+1,mStruct.nF);
            end
            TypeChecker.requireuniquechars({mStruct.bfaceData.id});
            This.bfaceData = mStruct.bfaceData;
         end
         
      catch ME
         if strcmp(ME.identifier,'MATLAB:nonExistentField')
            baseExc = MException('FVMLab:InsufficientData',...
               'The supplied data struct does not have all required fields.');
            ME = addCause(ME,baseExc);
         end
         rethrow(ME);
      end
        
    end % FvMesh constructor
    
                
end % PROTECTED METHODS



end
