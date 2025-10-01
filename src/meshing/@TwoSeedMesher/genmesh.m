function meshObj = genmesh(seedI,seedJ,boundarynames)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% FUNCTION meshObj = genmesh(seedI,seedJ)
%
% Creates a regular unstructured mesh from two non-colinear line seeds.
% Compared to the older genmeshV1 version, the vertices are also ordered as
% an internal block and a boundary vertex block (in accordance to the new
% specification). Further, the number scheme has changed and is now the
% following:
%     f22     f23     f24
% f15 c3  f3  c6  f6  c9  f18
%     f10     f11     f12
% f14 c2  f2  c5  f5  c8  f17       ^ J
%     f7      f8      f9            |
% f13 c1  f1  c4  f4  c7  f16       |
%     f19     f20     f21           o----> I
%
%  vertices follow the same pattern: first internal from S to N, and from W
%  to E, then boundary. The boundary vertices are numbered in the following
%  order: SW,NW,NE,SE, W, E, S, N
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


if nargin < 3 || isempty(boundarynames)
   boundarynames = {'WESTRAND','OOSTRAND','ZUIDRAND','NOORDRAND'};
end


% Error when seeds do not exactly share origin
if ~seedI.testShareOrigin(seedJ,0)
   % Tolerance is 0 --> very strict: vertex must be the same!!
   error('Supplied seeds must have exactly the same origin.')
end

% Error when I x J is not positive
crossp = seedI.displX*seedJ.displY - seedI.displY*seedJ.displX;
if crossp <= 0
   error('Supplied seeds must have positive cross product.')
end

% TODO: Error on colinearity of seeds.

nVertseedI = seedI.nVertices;
nVertseedJ = seedJ.nVertices;


nV = int32(nVertseedI*nVertseedJ);
nBv = int32(2*(nVertseedI-1) + 2*(nVertseedJ-1));
nIv = nV - nBv;

nF = seedI.nSegm*nVertseedJ + ...
   seedJ.nSegm*nVertseedI;
nBf = 2*(seedI.nSegm + seedJ.nSegm);
nIf = nF - nBf;

nPc = seedI.nSegm*seedJ.nSegm;

maxfaceconncells = int32(2);
maxfaceconnvertices = int32(2);


vCoord = zeros(2,nV);


faceFirstConnVertex = maxfaceconnvertices; % fixed connectivity: number instead of array
faceConnVertices = int32(zeros(1,maxfaceconnvertices*nF));

faceFirstConnCell = maxfaceconncells; % fixed connectivity: number instead of array
faceConnCells = int32(zeros(1,maxfaceconnvertices*nF));

cSeedI = [seedI.dX' ; seedI.dY'];
cSeedJ = [seedJ.dX' ; seedJ.dY'];
% seedIX = seedI.dX;
% seedIY = seedI.dY;
% seedJX = seedJ.dX;
% seedJY = seedJ.dY;

% Can be deleted:
originX = seedI.originX; % = seedJ.originX
originY = seedI.originY; % = seedJ.originY
vCoord(1,:) = originX; % vCoord(1,:) - 
vCoord(2,:) = originY; % vCoord(2,:) - 

% Generate boundary vertices and boundary faces============================
thisBVert = 0; % counters for creating boundary vertices
thisBFace = 0; % counters for creating boundary faces

% Create corner vertices: SW,NW,NE,SE
I = 1; J = 1; % SW
thisBVert = thisBVert + 1;
vCoord(1,nIv + thisBVert) = originX + cSeedI(1,I) + cSeedJ(1,J);      
vCoord(2,nIv + thisBVert) = originY + cSeedI(2,I) + cSeedJ(2,J);
   
I = 1; J = nVertseedJ; % NW
thisBVert = thisBVert + 1;
vCoord(1,nIv + thisBVert) = originX + cSeedI(1,I) + cSeedJ(1,J);      
vCoord(2,nIv + thisBVert) = originY + cSeedI(2,I) + cSeedJ(2,J);
   
I = nVertseedI; J = nVertseedJ; % NE
thisBVert = thisBVert + 1;
vCoord(1,nIv + thisBVert) = originX + cSeedI(1,I) + cSeedJ(1,J);      
vCoord(2,nIv + thisBVert) = originY + cSeedI(2,I) + cSeedJ(2,J);
   
I = nVertseedI; J = 1; % SE
thisBVert = thisBVert + 1;
vCoord(1,nIv + thisBVert) = originX + cSeedI(1,I) + cSeedJ(1,J);      
vCoord(2,nIv + thisBVert) = originY + cSeedI(2,I) + cSeedJ(2,J);

% Create W vertices I=1, J=2..nVJ-1 and faces
I = 1;
prevBVert = 1; % SW vertex
lastBVert = 2; % NW vertex
seedItmp(1) = cSeedI(1,I);
seedItmp(2) = cSeedI(2,I);
for J=2:nVertseedJ-1
   % Create new vertex
   thisBVert = thisBVert + 1;
   thisVert = nIv + thisBVert;
   vCoord(1,thisVert) = originX + seedItmp(1) + cSeedJ(1,J);
   vCoord(2,thisVert) = originY + seedItmp(2) + cSeedJ(2,J);
   % Create new face between previous vertex and the new vertex
   thisBFace = thisBFace + 1;
   faceConnVertices(2*(nIf + thisBFace - 1) + 1) = ...
      nIv + prevBVert;
   faceConnVertices(2*(nIf + thisBFace - 1) + 2) = ...
      nIv + thisBVert;
   prevBVert = thisBVert;
   % Set owner cell of this face
   ownerCell = (I+1-2)*(nVertseedJ-1) + (J - 1);
	faceConnCells(2*(nIf + thisBFace - 1) + 1) = ownerCell;
end
% Create new face between previous vertex and the last vertex
thisBFace = thisBFace + 1;
faceConnVertices(2*(nIf + thisBFace - 1) + 1) = ...
	nIv + prevBVert;
faceConnVertices(2*(nIf + thisBFace - 1) + 2) = ...
	nIv + lastBVert;
% Set owner cell of this face
ownerCell = (I+1-2)*(nVertseedJ-1) + (nVertseedJ - 1);
faceConnCells(2*(nIf + thisBFace - 1) + 1) = ownerCell;


% Create E vertices I=nVI, J=2..nVJ-1 and faces
I = nVertseedI;
prevBVert = 4; % SE vertex
lastBVert = 3; % NE vertex
seedItmp(1) = cSeedI(1,I);
seedItmp(2) = cSeedI(2,I);
for J=2:nVertseedJ-1
   % Create new vertex
   thisBVert = thisBVert + 1;
   thisVert = nIv + thisBVert;
   vCoord(1,thisVert) = originX + seedItmp(1) + cSeedJ(1,J);
   vCoord(2,thisVert) = originY + seedItmp(2) + cSeedJ(2,J);
   % Create new face between previous vertex and the new vertex
   thisBFace = thisBFace + 1;
   faceConnVertices(2*(nIf + thisBFace - 1) + 1) = ...
      nIv + thisBVert;
   faceConnVertices(2*(nIf + thisBFace - 1) + 2) = ...
      nIv + prevBVert;      
   prevBVert = thisBVert;
   % Set owner cell of this face
   ownerCell = (I-2)*(nVertseedJ-1) + (J - 1);
	faceConnCells(2*(nIf + thisBFace - 1) + 1) = ownerCell;   
end
% Create new face between previous vertex and the last vertex
thisBFace = thisBFace + 1;
faceConnVertices(2*(nIf + thisBFace - 1) + 1) = ...	
   nIv + lastBVert;
faceConnVertices(2*(nIf + thisBFace - 1) + 2) = ...
   nIv + prevBVert;
% Set owner cell of this face
ownerCell = (I-2)*(nVertseedJ-1) + (nVertseedJ - 1);
faceConnCells(2*(nIf + thisBFace - 1) + 1) = ownerCell;

% Create S vertices I=2..nVI-1, J=1 and faces
J = 1;
prevBVert = 1; % SW vertex
lastBVert = 4; % SE vertex
seedJtmp(1) = cSeedJ(1,J);
seedJtmp(2) = cSeedJ(2,J);
for I=2:nVertseedI-1
	% Create new vertex
   thisBVert = thisBVert + 1;
   thisVert = nIv + thisBVert;
   vCoord(1,thisVert) = originX + cSeedI(1,I) + seedJtmp(1);
   vCoord(2,thisVert) = originY + cSeedI(2,I) + seedJtmp(2);
   % Create new face between previous vertex and the new vertex
   thisBFace = thisBFace + 1;
   faceConnVertices(2*(nIf + thisBFace - 1) + 1) = ...
      nIv + thisBVert;
   faceConnVertices(2*(nIf + thisBFace - 1) + 2) = ...
      nIv + prevBVert;      
   prevBVert = thisBVert;
   % Set owner cell of this face   
	ownerCell = (I-2)*(nVertseedJ-1) + (J + 1) - 1;
   faceConnCells(2*(nIf + thisBFace - 1) + 1) = ownerCell;
end
thisBFace = thisBFace + 1;
faceConnVertices(2*(nIf + thisBFace - 1) + 1) = ...
	nIv + lastBVert;
faceConnVertices(2*(nIf + thisBFace - 1) + 2) = ...
   nIv + prevBVert;
% Set owner cell of this face   
ownerCell = (nVertseedI-2)*(nVertseedJ-1) + (J + 1) - 1;
faceConnCells(2*(nIf + thisBFace - 1) + 1) = ownerCell;


% Create N vertices I=2..nVI-1, J=nVJ and faces
J = nVertseedJ;
prevBVert = 2; % NW vertex
lastBVert = 3; % NE vertex
seedJtmp(1) = cSeedJ(1,J);
seedJtmp(2) = cSeedJ(2,J);
for I=2:nVertseedI-1
	% Create new vertex
   thisBVert = thisBVert + 1;
   thisVert = nIv + thisBVert;
   vCoord(1,thisVert) = originX + cSeedI(1,I) + seedJtmp(1);
   vCoord(2,thisVert) = originY + cSeedI(2,I) + seedJtmp(2);
   % Create new face between previous vertex and the new vertex
   thisBFace = thisBFace + 1;
   faceConnVertices(2*(nIf + thisBFace - 1) + 1) = ...
      nIv + prevBVert;
   faceConnVertices(2*(nIf + thisBFace - 1) + 2) = ...
      nIv + thisBVert;
   prevBVert = thisBVert;
   % Set owner cell of this face   
	ownerCell = (I-2)*(nVertseedJ-1) + nVertseedJ - 1;
   faceConnCells(2*(nIf + thisBFace - 1) + 1) = ownerCell;
end
thisBFace = thisBFace + 1;
faceConnVertices(2*(nIf + thisBFace - 1) + 1) = ...
	nIv + prevBVert;
faceConnVertices(2*(nIf + thisBFace - 1) + 2) = ...
	nIv + lastBVert;
% Set owner cell of this face   
ownerCell = (nVertseedI-2)*(nVertseedJ-1) + nVertseedJ - 1;
faceConnCells(2*(nIf + thisBFace - 1) + 1) = ownerCell;


%Generate internal vertices================================================
thisVert = 0;
for I=2:nVertseedI-1
   for J=2:nVertseedJ-1
      thisVert = thisVert + 1;
      vCoord(1,thisVert) = originX + cSeedI(1,I) + cSeedJ(1,J);
      vCoord(2,thisVert) = originY + cSeedI(2,I) + cSeedJ(2,J);
   end
end

% Generate internal faces and the cells====================================
thisFace = 0;
% (1) Generate internal vertical faces and link to cell
thisVert = 0;
for I=2:nVertseedI-1
   prevVert = nIv + 4 + 2*(nVertseedJ-2) + (I-1);
   lastVert = prevVert + (nVertseedI-2);
   for J=2:nVertseedJ-1
      thisVert = thisVert + 1;
      % Create new face between previous vertex and the new vertex
      thisFace = thisFace + 1;
      faceConnVertices(2*(thisFace - 1) + 1) = ...
         thisVert;
      faceConnVertices(2*(thisFace - 1) + 2) = ...
         prevVert;
      prevVert = thisVert;
      % Set owner and neighbor cells of this face
      ownerCell = (I-2)*(nVertseedJ-1) + J - 1;
      faceConnCells(2*(thisFace - 1) + 1) = ownerCell;
      faceConnCells(2*(thisFace - 1) + 2) = ownerCell + (nVertseedJ-1);
   end
   % Create new face between previous vertex and the last vertex
   thisFace = thisFace + 1;
   faceConnVertices(2*(thisFace - 1) + 1) = ...
      lastVert;
   faceConnVertices(2*(thisFace - 1) + 2) = ...
      prevVert;
	% Set owner and neighbor cells of this face
   ownerCell = (I-2)*(nVertseedJ-1) + nVertseedJ - 1;
   faceConnCells(2*(thisFace - 1) + 1) = ownerCell;
   faceConnCells(2*(thisFace - 1) + 2) = ownerCell + (nVertseedJ-1);   
end


% (2) Generate internal horizontal faces and link to cell
for J=2:nVertseedJ-1
   prevVert = nIv + 4 + (J-1);
   lastVert = prevVert + (nVertseedJ-2);
   thisVert = 1 - (nVertseedJ-2) + (J-2);
   for I=2:nVertseedI-1
      %thisVert = thisVert + (nVertseedI-2);
      thisVert = thisVert + (nVertseedJ-2);
      % Create new face between previous vertex and the new vertex
      thisFace = thisFace + 1;
      faceConnVertices(2*(thisFace - 1) + 1) = ...
         prevVert;
      faceConnVertices(2*(thisFace - 1) + 2) = ...
         thisVert;
      prevVert = thisVert;
      % Set owner and neighbor cells of this face   
      ownerCell = (I-2)*(nVertseedJ-1) + J - 1;
      faceConnCells(2*(thisFace - 1) + 1) = ownerCell;
      faceConnCells(2*(thisFace - 1) + 2) = ownerCell + 1;      
   end
   % Create new face between previous vertex and the last vertex
   thisFace = thisFace + 1;
   faceConnVertices(2*(thisFace - 1) + 1) = ...
      prevVert;
   faceConnVertices(2*(thisFace - 1) + 2) = ...
      lastVert;
	% Set owner and neighbor cells of this face   
   ownerCell = (nVertseedI-2)*(nVertseedJ-1) + J - 1;
   faceConnCells(2*(thisFace - 1) + 1) = ownerCell;
   faceConnCells(2*(thisFace - 1) + 2) = ownerCell + 1;
end

% TODO: move this code to domain creation
% % Construct cellFirstConnFace & cellConnFaces==============================
% cellnFaces = int32(4);
% cellFirstConnFace = (1 + cellnFaces*( (1:1:nPc+1) - 1));
% cellConnFaces = int32(zeros(1,cellnFaces*nPc));
% for jFace = 1:nF
%    ownerCell = faceConnCells(faceFirstConnCell*(jFace-1) + 1);
%    basePtr = cellnFaces*(ownerCell-1) + 1;
% %    fprintf('%g, %g\n',jFace,ownerCell)
%    while cellConnFaces(basePtr) ~= 0
%       basePtr = basePtr + 1;
%    end
%    cellConnFaces(basePtr) = jFace;
%    
%    neighborCell = faceConnCells(faceFirstConnCell*(jFace-1) + 2);
%    if neighborCell > 0
%       basePtr = cellnFaces*(neighborCell-1) + 1;   
%       while cellConnFaces(basePtr) ~= 0
%          basePtr = basePtr + 1;
%       end
%       cellConnFaces(basePtr) = jFace;
%    end   
% end


% Create data structure for boundary faces ========================
bFaceData = struct('id',{},'range',{});

endFace = nIf;

iBoundary = 1;
startFace = endFace + 1;
endFace = startFace + seedJ.nSegm - 1;
bFaceData(iBoundary).id = boundarynames{iBoundary};
bFaceData(iBoundary).range = [startFace ; endFace];

iBoundary = iBoundary + 1;
startFace = endFace + 1;
endFace = startFace + seedJ.nSegm - 1;
bFaceData(iBoundary).id = boundarynames{iBoundary};
bFaceData(iBoundary).range = [startFace ; endFace];

iBoundary = iBoundary + 1;
startFace = endFace + 1;
endFace = startFace + seedI.nSegm - 1;
bFaceData(iBoundary).id = boundarynames{iBoundary};
bFaceData(iBoundary).range = [startFace ; endFace];

iBoundary = iBoundary + 1;
startFace = endFace + 1;
endFace = startFace + seedI.nSegm - 1;
bFaceData(iBoundary).id = boundarynames{iBoundary};
bFaceData(iBoundary).range = [startFace ; endFace];

% Put everything in a struct ==================================
mStruct.dim    = int32(2);
mStruct.nV     = nV;
mStruct.nF     = nF;
mStruct.nPc    = nPc;
mStruct.nIv    = nIv;
mStruct.nBv    = nBv;
mStruct.nIf    = nIf;
mStruct.nBf    = nBf;

mStruct.vCoord       = vCoord;

mStruct.fNbVLoc      = faceFirstConnVertex;
mStruct.fNbV         = faceConnVertices;
mStruct.fNbCLoc      = faceFirstConnCell;
mStruct.fNbC         = faceConnCells;
mStruct.maxNoCNbF    = int32(4);
mStruct.maxNoCNbV    = int32(4);

% This data will be calculated on domain creation
% mStruct.cellToFacePtrFirst = cellFirstConnFace;
% mStruct.cellToFaceList = cellConnFaces;

mStruct.noNamedBfaceZones ...
                     = cast(length(bFaceData),'int32');
mStruct.bfaceData    = bFaceData;

% Create mesh from minimal data ===============================
% minimal data is face definition & cell definition
meshObj = FvMesh.fromstruct(mStruct);

end