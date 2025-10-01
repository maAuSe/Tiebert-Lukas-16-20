classdef FvDomain < FvmLabEntity
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% CLASS FvDomain 
%
% A domain is the basic entity on which field problems are solved. It holds
% all raw topological and geometrical data required for carrying out finite
% volume discretisation on primal cells. Further a set of higher level
% objects is maintained by the domain in order to assist in a higher
% level coding approach. These include default vertex, face and cell zones,
% together with (mesh- or user-) defined boundary zones.
%
% Construction of a FvDomain is based on a FvMesh which holds only a
% minimal data set for describing the topology and geometry. Optionally, an
% id can be supplied for the domain:
%     FvDomainObj = FvDomain(MeshObj,domName)
% 
% NOTE: 
%   Every boundary face has a ghost cell connected at the outside of the
%   physical domain. These ghost cells are not *real* cells in the sense
%   that they do not result from a closed set of faces. They only have a
%   centroid (node location) and a volume. The volume is the same as the n1
%   cell of the boundary face. The ghostcell centroid is obtained as the
%   n1 centroid point-mirrorred through the boundary face centroid.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



properties (GetAccess=public,SetAccess=protected)

   % IMPORTANT VALUES -----------------------------------------------------
   dim         % int32
   % Spatial dimension of mesh
   
   
   % TOPOLOGICAL DATA -----------------------------------------------------   
   nV          % int32
   % Number of Vertices  
   nIv         % int32
   % Number of Internal vertices
   nBv         % int32
   % Number of Boundary vertices
   
   nF          % int32
   % Number of Faces
   nIf         % int32
   % Number of Internal faces
   nBf         % int32
   % Number of Boundary faces
   
   nPc         % int32
   % Number of dom (physical) Cells
	nGc      % int32
   % Number of Ghost cells
   nC       % int32
   % Total number of cells nC = nPc + nGc
      
   fNbVLoc     % int32 == 2
   fNbV        % int32[fNbVLoc*nF]
   % Face neighbor vertices
   
   fNbCLoc     % int32 == 2   % Always (1D, 2D, 3D)
   fNbC        % int32[fNbCLoc*nF]
   % Face neighbor cells.
   % Internal faces connect to c1=pC1 and c2=pC2.
   % Boundary faces connect to c1=pC1 and c2=gC2+nPc.
   
   maxnoCNbF   % int32
   maxnoCNbV   % int32
   % Maximum number of connections in the mesh. Can be used to assist
   % algorithms in allocating sufficient memory and setting up data
   % structures. If it is 0, the maximum connectivity is unknown.
   
   vNbFLoc     % int32[nV+1]
   vNbF        % int32[2*nF]
   % Vertex neighbor faces. Not ordered
   
   cNbFLoc     % int32[nC+1]
   cNbF        % int32[2*nIf+nBf]
   % Cell neighbor faces. Ordered (counter-)clock-wise
   
   vNbCLoc     % int32[nV+1]
   vNbC        % int32[???]   % 2*nIf+nBf
   % Vertex neighbor cells. Not ordered
   
   cNbVLoc     % int32[nC+1]
   cNbV        % int32[2*nIf + nBf]
   % Cell neighbor vertices. Ordered (counter-)clock-wise for physical
   % cells.
   
   
   % GEOMETRICAL DATA -----------------------------------------------------
   vCoord     	% double[2,nV]
   %  Vertices (X,Y)-coordinates
   fCoord      % double[2,nF]
   %  Faces centroid (X,Y)-coordinates
   cCoord      % double[2,nC]
   %  Cells centroid (X,Y)-coordinates
   
   fArea       % double[1,nF]
   %  Faces areas  
   fNormal     % double[2,nF]
   %  Faces unit normal (X,Y)-coordinates
   %  The normal vector points outside the n1 cell and inside the n2 cell:
   %     (1) dot(fNormal,fTangent) == 0;
   %  	(2) dot(fNormal,fXi) > 0
   %     (3) norm(fNormal) == 1
   fTangent    % double[2,nF]   
   %  Faces unit tangent (X,Y)-coordinates
   %  The tangent vector is defined as follows:
   %     (1) cross(fTangent,fNormal) is in positive z-direction.
   %     (2) norm(fTangent) == 1      
         
   cVol        % double[1,nC]
   %  Cells volumina   
   fXi         % double[2,nF]
   %  vector pointing from n1 cCoord to n2 cCoord   
   fXiMag      % double[2,nF]
   % Magnitude of Xi vector
   fXiLambda   % double[2,nF]
   % fraction of fXi until Xi intersects the face.
   % Measured, starting from N1
   fDs         % double[2,2*nF]
   %  Is indexed from fNbCLoc == 2.
   %  (1) The first item is vector pointing from n1 cCoord to fCoord.
   %  (2) The second item is vector pointing from n2 cCoord to fCoord.
   
         
   % HIGH LEVEL ZONE DATA -------------------------------------------------
   allCells       % Default zone containing all cells
   allFaces       % Default zone containing all faces
   allBfaces      % Default zone containing all boundary faces
   allVerts       % Default zone containing all vertices
   
   ZoneList       % List container
   % A list container with unique zone objects. Zones can be queried by their
   % unique "id".
   

end % PROPERTIES




methods (Access=public)

    
	function this = FvDomain(MeshObj,domName)
      
      this = this@FvmLabEntity(domName);
      
      TypeChecker.requireclass(MeshObj,'FvMesh');
      
      % Extract raw data from mesh ========================================
      this.dim          = MeshObj.dim;
      
      this.nV           = MeshObj.nV;
      this.nIv          = MeshObj.nIv;
      this.nBv          = MeshObj.nBv;      
      
      this.nF           = MeshObj.nF;
      this.nIf          = MeshObj.nIf;
      this.nBf          = MeshObj.nBf;

      this.nPc          = MeshObj.nPc;
      this.nGc          = this.nBf; % Note: depends on type of gCells
      this.nC           = this.nPc + this.nGc;
      
      this.vCoord       = MeshObj.vCoord;
      
      this.fNbVLoc      = MeshObj.fNbVLoc;
      this.fNbV         = MeshObj.fNbV;
      this.fNbCLoc      = MeshObj.fNbCLoc;
      this.fNbC         = MeshObj.fNbC;
      this.maxnoCNbF    = MeshObj.maxNoCNbF;
      this.maxnoCNbV    = MeshObj.maxNoCNbV;
      
      % Add mirror ghostcells to fNbC:
      this.fNbC = add_mirrgc_to_fNbC(...
         this.fNbCLoc,this.fNbC,this.nIf,this.nBf,this.nPc);
                  
      % Calculate other topological data ==================================      
      [this.vNbFLoc this.vNbF] = calc_vNbF(...
         this.nV,this.nF,this.fNbVLoc,this.fNbV);
      
      [this.cNbFLoc this.cNbF] = calc_cNbF(...
         this.nC,this.nPc,this.nF,this.maxnoCNbF,...
         this.fNbCLoc,this.fNbC,this.fNbVLoc,this.fNbV);
      
      [this.cNbVLoc this.cNbV] = calc_cNbV(...
         this.nC,this.nPc,this.maxnoCNbV,...
         this.cNbFLoc,this.cNbF,this.fNbVLoc,this.fNbV);
                  
      [this.vNbCLoc this.vNbC] = inverttopology(...
         this.nC,this.nV,this.cNbVLoc,this.cNbV); % TODO: validate this routine 


      % Calculate other geometrical data ==================================
      this.fCoord = calc_fCentr(...
         this.nF,this.fNbVLoc,this.fNbV,this.vCoord);

      this.fArea = calc_fArea(this.nF,...
         this.fNbVLoc,this.fNbV,this.vCoord);
    
      this.fNormal = calc_fNormal(this.fNbVLoc,this.fNbV,...
         this.vCoord,this.fArea);
      
      this.fTangent = calc_fTangent(this.nF,this.fNormal);
             
      [this.cVol this.cCoord] = calc_cVol_cCoord(this.nC,this.nPc,...
         this.vCoord,this.fCoord,...
         this.cNbFLoc,this.cNbF,this.fNbVLoc,this.fNbV,...
         this.fNbCLoc,this.fNbC);
           
      this.fXi = calc_fXi(this.nF,...
         this.cCoord,...
         this.fNbCLoc,this.fNbC);
      this.fXiMag = calc_fXiMag(this.nF,this.fXi);
      this.fXiLambda = calc_fXiLambda(this.nF,...
         this.fXi,this.fTangent,this.cCoord,this.vCoord,...
         this.fNbCLoc,this.fNbC,this.fNbVLoc,this.fNbV);
      
      this.fDs = calc_fDs(this.nF,...
         this.fCoord,this.cCoord,...
         this.fNbCLoc,this.fNbC);
      
      
      % Create higher-level objects =======================================
      ZoneList = FvmLabEntityList(); % Entitylist without duplicate id's.
      this.ZoneList = ZoneList;
      
      % Construct default zones      
      this.allVerts  = newzone('allVerts', this,[1;this.nV], 'Vert');
      this.allFaces  = newzone('allFaces', this,[1;this.nF], 'Face');
      this.allBfaces = newzone('allBfaces',this,[this.nIf+1;this.nF],'Bface');
      this.allCells  = newzone('allCells', this,[1;this.nC], 'Cell');
      
      
      add(ZoneList,this.allVerts);
      add(ZoneList,this.allFaces);
      add(ZoneList,this.allBfaces);
      add(ZoneList,this.allCells);      
      
      % Construct named boundary face zones as defined in the mesh
      bfaceData = MeshObj.bfaceData;
      for jBzone = 1:MeshObj.noNamedBfaceZones         
         add(ZoneList,...
             newzone(bfaceData(jBzone).id,this,...
                     bfaceData(jBzone).range,'Bface'));                                        
      end

   end % FvDomain constructor
   
   
   display(this)
   ZoneObj = getzone(this,id)
   undefinezone(this,id)

    
end % METHODS




end



