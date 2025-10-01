function require_validrange(range,varargin)

   try
   	TypeChecker.requireclass(range,'int32');
      TypeChecker.requirenrows(range,2);        
	catch ME
      throwAsCaller(ME);
	end

   for iCol = 1:size(range,2)
      diff = range(2,iCol) - range(1,iCol);
      if diff < 0 || ...
            (iCol > 1 && range(1,iCol) <= range(2,iCol-1))
         ME = MException('FVMLab:TypeChecker:RangeViolation',...
            ['\nFVMLab:TypeChecker:RangeViolation',...
             '\nRange arrays must be monotonic increasing'...
             '\n']);
         throwAsCaller(ME);
      end
   end
   
   if nargin >= 2 && min(min(range)) < varargin{1}
      ME = MException('FVMLab:TypeChecker:RangeViolation',...
         ['\nFVMLab:TypeChecker:RangeViolation',...
          '\nThe supplied range data defines a range lower than the allowed minimum of %g'...
          '\n'],varargin{1});
      throwAsCaller(ME);
   end
   
   if nargin >= 3 && max(max(range)) > varargin{2}
      ME = MException('FVMLab:TypeChecker:RangeViolation',...
         ['\nFVMLab:TypeChecker:RangeViolation',...
          '\nThe supplied range data defines a range higher than the allowed maximum of %g'...
          '\n'],varargin{2});
      throwAsCaller(ME);
   end
    
    
    
end