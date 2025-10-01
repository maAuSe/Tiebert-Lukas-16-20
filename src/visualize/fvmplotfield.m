function gh = fvmplotfield(a,scaleSpec,lw,c1,c2)
%FVMPLOTFIELD plots a cell field or a vertex field
%
%  FVMPLOTFIELD(MYFIELD,SCALE,LW,C1,C2)
%     MYFIELD is a vertex field or cell field
%     SCALE is 'lin' or 'log'
%     LW is the line width in pts of the mesh lines. When lw==0, no mesh is plotted
%     C1 specifies the first component to be plotted. Only required for
%     vector fiels or tensor fields.
%     C2 specifies the second component to be plotted. Only required for 
%     tensor fields.
%
%  NOTE: In case of a cell Field, the ghost cell values are not plotted.
%


   % lw is LineWdth in points
   lw = abs(lw);
   if lw == 0
      ec = 'none';
      ls = 'none';
      lw = 1;
   else
      ls = '-';
      ec = [1 1 1];
   end
   
   
   if ~isa(a,'Field')
      error('A Field must be supplied in order to create a field plot.')
   end

   if ~haselements(a)
      error('The supplied field is empty. Cannot make the field plot.')
   end
   
   % Retrieve the desired component for plotting
   if isscalar(a)
      fval = get(a);
      fval = fval'; % Transpose is for interfacing with Matlab plot routines.
   elseif isvector(a)
      fval = get(a,c1);
      fval = fval'; % Transpose is for interfacing with Matlab plot routines.
   else % istensor
      fval = get(a,c1,c2);
      fval = fval'; % Transpose is for interfacing with Matlab plot routines.
   end

   
   zone = getzone(a);
   dom  = getdomain(a);
   
   % Retrieve Mesh vertex coordinate values   
   vcoords = dom.vCoord;
   vcoords = vcoords'; % Transpose is for interfacing with Matlab plot routines.
   
   % Calculate connection matrix for patch faces definition
   % NOTE: patch and face are Matlab plot routines terminology.
   nPc   = dom.nPc;
   maxnoCNbV = dom.maxnoCNbV;
   cNbVLoc   = dom.cNbVLoc;
   cNbV      = dom.cNbV;
   faceDef = zeros(nPc,maxnoCNbV); % initialize   
   nextPtr = 1;
   for jC = 1:nPc
      currPtr = nextPtr;
      nextPtr = cNbVLoc(jC + 1);
      thisCNoNbV = nextPtr - currPtr;      
      for k = 1:thisCNoNbV
         faceDef(jC,k) = cNbV(currPtr+k-1);
      end
      % PRE: There must be at least three vertex neighbors!!
      for k = thisCNoNbV+1:maxnoCNbV
         faceDef(jC,k) = faceDef(jC,k-1);
      end
   end
   
   
   if isa(zone,'CellZone')
      % Cells are given the color corresponding to the cell values
      % The values in fval are ordered as the cells.      
      fc = 'flat';
      fval = fval(1:nPc);
   elseif isa(zone,'VertZone')
      % Cells are given an interpolated color corresponding to the vertex
      % values. The values in fval are ordered as the vertices. ?faces?.      
      fc = 'interp';   
   else
      error('Only CellZone or VertZone fields can be plotted')
   end
      

   if (nargin > 1) && strcmp(scaleSpec,'log')
      fval = log(abs(fval))/log(10);
   end
   
   
   % Create patch and set properties
%    gh = patch('Vertices',vcoords,'Faces',faceDef,'FaceVertexCData',fval,...
%          'NormalMode','manual','HitTest','on',...
%          'FaceColor',fc,...
%          'LineStyle',ls,'LineWidth',lw,'EdgeColor',ec);

    gh = patch('Vertices',vcoords,'Faces',faceDef,'FaceVertexCData',fval,...
         'VertexNormalsMode','manual','FaceNormalsMode','manual',...
         'HitTest','on','FaceColor',fc,...
         'LineStyle',ls,'LineWidth',lw,'EdgeColor',ec);
   
    
end