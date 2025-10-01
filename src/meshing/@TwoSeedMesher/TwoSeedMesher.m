classdef TwoSeedMesher < handle
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% CLASS TwoSeedMesher
%
% This class provides a static method to generate of mesh from two line
% seeds which share the same origin.
%
%   meshObj = genmesh(seedI,seedJ)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    

properties (Constant,GetAccess=private)
   
   TOL = 1e-14;
   
end % PRIVATE CONSTANTS




methods (Static=true)
   
   meshObj = genmeshV1(seedI,seedJ)             % TODO: Add optional domain name
   meshObj = genmesh(seedI,seedJ,boundarynames) % TODO: Add optional domain name
   
end % STATIC METHODS




end % TwoSeedMesher