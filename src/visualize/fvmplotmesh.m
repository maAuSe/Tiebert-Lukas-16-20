function gh = fvmplotmesh(dom,lw)
%FVMPLOTMESH plots the mesh of a domain.
%
%  FVMPLOTMESH(DOM,LW)
%     DOM is a FvDomain
%     LW is the line width in pts of the mesh lines.
%
%  NOTE: When fields are plotted, there is no need to call this routine 
%  first in order to plot the mesh. Calling fvmplotfield with LW>0 plots
%  both field and mesh.
%

   DummyField = Field(dom.allVerts,0);
   reset(DummyField,NaN);
  
   gh = fvmplotfield(DummyField,'lin',lw);

end