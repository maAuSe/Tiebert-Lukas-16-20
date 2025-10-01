classdef ScalarFvEqn2 < handle
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% CLASS ScalarFvEqn2
%
% Class that represents a system of linear equations A*x=b in an
% intermediate form that is determined by the topology of a FvDomain.
%
% The n by n matrix A is represented by the adata vector of size nnz
% holding all possible coefficients of a compact discretization. The
% following identities apply:
% A(i,i) == adata(i)
%     for every cell: i in [1..dom.nC]
% A(fn1(j),fn2(j)) == adata(n+2*j-1);
%     for every internal face: j in [1..dom.nIf]
% A(fn2(j),fn1(j)) == adata(n+2*j);
%     for every internal face: j in [1..dom.nIf]
% A(fn1(j),fn2(j)) == adata(n+2*dom.nIf+2*j-1)
%     for every boundary face: j in [dom.nIf+1..dom.nF]
% A(fn2(j),fn1(j)) == adata(n+2*dom.nIf+2*j)
%     for every boundary face: j in [dom.nIf+1..dom.nF]
% Where fn1(j) and fn2(j) are the first and second neighbor cells of face j.
% Note that for boundary faces fn2(j) refers to a ghostcell which is used
% for imposing boundary conditions.
% Thus adata holds:
% First the diagonal elements of A;
% Then, off-diagonal elements related to internal face contributions in 
%       fn1 and fn2 equations;
% Then, off-diagonal elements related to boundary face contributions in 
%       fn1 and fn2 equations;
%
% The right-hand-side of the system b is equal to the dense vector bdata
%
% A ScalarFvEqn can be converted to a matlab sparse coefficient matrix and
% a dense right-hand-side in order to interface the matlab built-in linear
% solvers. In order to do this, one must set the adata and bdata properties
% and call [A,b] = to_msparse(ScalarFvEqnObj);
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




properties (GetAccess=public,SetAccess=public)

	adata         % double mArray(nnz,1)
   % Coefficients in matrix A.
    
   bdata         % double mArray(n,1)
   % Right-hand-side b
   
end


properties (GetAccess=public,SetAccess=private)
    
   n              % double mArray(1,1) 
   % Problem size
   
   nnz            % double mArray(1,1) 
   % Number of nonzeros:
   % nnz = nPCells + nGCells + 2*nIfaces + 2*nBfaces      
     
   ir            % double mArray(nnz,1)
   jc            % double mArray(nnz,1)
   % Row and culumn indices of the adata elements.
   % Forming the msparse is then easy: A = sparse(ir,jc,adata,n,n)
       
   dom           % FvDomain
    
    
end



methods (Access=public)
    
    function this = ScalarFvEqn2(dom)
    % CONSTRUCTOR for ScalarFvEqn    
    % RETURNS: a ready-to-use ScalarFvEqn object holding the trivial 
    % equation 0*x==0
    
      TypeChecker.requireclass(dom,'FvDomain')
      this.dom = dom;
        
      nC    = dom.nC;
      nIf   = dom.nIf;
      nBf   = dom.nBf;
      
      n     = nC;
      nnz   = nC + 2*nIf + nBf + nBf;
            
      this.n   = n;
      this.nnz = nnz;

      fNbC     = dom.fNbC;
      fNbCLoc  = dom.fNbCLoc;
      if fNbCLoc ~= 2
         error('The face neighbor cells are not conform the specification')
      end
      
      [this.ir this.jc] = ScalarFvEqn2.calc_ir_jc(nnz,nIf+nBf,fNbC,n);
        
      reset(this);            
        
    end
    
    reset(this)
    % Reset to trivial equation 0*x==0
    
    [A,b] = to_msparse(this)
    % Convert to matlab sparse linear system
    
    
end % METHODS




methods (Static=true,Access=private)
    
   [ir,jc] = calc_ir_jc(nnz,nF,nIf,fNbC,n);
    % Computes the row and column indices of adata elements.    
    
end % STATIC methods
    
    
end