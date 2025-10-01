classdef LineSeed < handle
   
   %----------------------------------------------------------------------
   %  Private fields
   %----------------------------------------------------------------------
   properties (SetAccess=private)
      
      nVertices   % int32(1,1)
      %   Number of vertices in seed
      nSegm       % int32(1,1)
      %   Number of line segments in seed
      dX           % double(nVertices,1)
      %   X-coordinates of seed vertices relative to seed origin
      dY           % double(nVertices,1)
      %   Y-coordinates of seed vertices relative to seed origin
      originX     % double(1,1)
      %   X-coordinate of origin vector (0 -> start)
      originY     % double(1,1)
      %   Y-coordinate of origin vector (0 -> start)
      displX      % double(1,1)
      %   X-coordinate of displacement vector (start -> end)
      displY      % double(1,1)
      %   Y-coordinate of displacement vector (start -> end)
      
   end % properties
   
   
   
   %----------------------------------------------------------------------
   %  Public methods
   %----------------------------------------------------------------------
   methods
      
      function obj = LineSeed()
         
         obj.nVertices = int32(0);
         obj.nSegm = int32(0);
         obj.dX = zeros(obj.nVertices,1);
         obj.dY = zeros(obj.nVertices,1);
         obj.originX = 0;
         obj.originY = 0;
         obj.displX = 0;
         obj.displY = 0;
         
         
      end % LineSeed
      
      
      function [oX,oY] = getOrigin(obj)
         
         oX = obj.originX;
         oY = obj.originY;
         
      end % getOrigin
      
      
      function [dX,dY] = getDestination(obj)
         
         dX = obj.dX(end);
         dY = obj.dY(end);
         
      end % getDestination
      
      
      function test = testShareOrigin(obj,lineSeed,tol)
         
         [oX oY] = obj.getOrigin();
         [lX lY] = lineSeed.getOrigin();
         
         dist = sqrt((oX - lX)^2 + (oY - lY)^2);
         
         if dist > tol
            test = false;
         else
            test = true;
         end
         
      end
      
      
      function test = isOrthogonalTo(obj,lineSeed,tol)
         
         innerProd = obj.displX * lineSeed.displX + ...
            obj.displY * lineSeed.displY;
         
         if abs(innerProd) > tol
            test = false;
         else
            test = true;
         end
         
         
      end
      
      
      
   end % Public methods --------------------------------------------------
   
   
   
   %----------------------------------------------------------------------
   %  Static public methods
   %----------------------------------------------------------------------
   methods (Static=true,Access=public)
      
      
      obj = lineSeedOneWayBias(originV,displV,...
            nParts,expFac,biasStart)
               
      
   end % Static public methods -------------------------------------------
   
   
   
   %----------------------------------------------------------------------
   %  Static protected methods
   %----------------------------------------------------------------------
   methods (Static=true,Access=protected)
      
      
      function length = expansionLength(expFac,nParts)
         
         % %number must be > 1
         % if number < 1
         %     error('number must be >= 1')
         % end
         
         length = 0;
         for i = 1:1:nParts
            length = length + (expFac)^(i-1);
         end
         
      end %expansionLength
      
   end % Static protected methods ----------------------------------------
   
   
end % LineSeed