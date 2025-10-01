%==========================================================================
%
% Example solver using the FVMLab framework 4 students
%
% Purpose: Provides code structure for solving a scalar conservation
%          equation using data structures provided by the framework.
%
% by Frederik Rogiers
%
%==========================================================================
function result = examplesolver(casedef)

dom = casedef.dom;

% Create field objects
T = Field(dom.allCells,0);      % Temperature [K] (scalar); empty field
reset(T,0);                     % Reset with all zeros

% ... Create all other required data structures

% Create an equation object for holding a scalar conservation equation
eqn = ScalarFvEqn2(dom);

iterate = true;
niter = 0;
while iterate   
   
   niter = niter+1;
   
   % Set all terms to zero
   reset(eqn); 
   
   % Compute coefficients for physical cell eqns and add them to eqn object
   % ...
   % Compute coefficients for ghost cell eqns and add them to eqn object
   % ...
   eqn.adata = rand(eqn.nnz,1); % just a meaningless example adata
   eqn.bdata = rand(eqn.n,1);   % just a meaningless example bdata
      
   % Create a matlab sparse linear system from the eqn object
   [A,b] = to_msparse(eqn);
   x = get(T);
   x = x';
   
   % Check tolerance and iteration count
   TRes = b-A*x;
   TResnorm = norm(TRes);         
   if TResnorm < casedef.iteration.TTol
      Tconverged = true;
      iterate = false;
   elseif niter > casedef.iteration.maxniter
      Tconverged = false;
      iterate = false;
%    elseif checkstoprequest(stopmon)
%       Tconverged = false;
%       iterate = false;
   else
      x = A\b; % Direct sparse solver.
               % Alternatives: gmres, bicgstabb, ...
      set(T,x'); % Put algebraic solution in the Field
   end
      
   
end % iterate

result.endtime = now; % call datestr(now) for displaying this time 
result.Tconverged = Tconverged;
result.niter = niter;
result.TResnorm = TResnorm;
result.TRes = Field(dom.allCells,0);
   set(result.TRes,TRes');
result.T = T;


end



