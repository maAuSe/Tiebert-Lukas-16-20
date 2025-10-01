classdef IZone < handle
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% INTERFACE IZOne
%
% This interface desribes the ability of subtypes to return the number of
% elements the zone must contain, their coordinates, and the domain the
% zone belongs to.
% Data defined on a zone thus refers to a set of discrete locations defined
% by the subtype.
% 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




methods (Access=public,Abstract=true)
       
   DomObj = getdomain(This)      % returns the domain
   dim = getdim(This)            % returns the dimension of the domain
   n = getelcount(This)          % returns the number of elements in the zone
   cdata = getcoords(This)       % returns the coordinates of elements
   cdata = getallcoords(This)    % returns all coordinates of alike elements 
   
end % ABSTRACT METHODS




end










