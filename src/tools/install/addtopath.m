function addtopath(somePath,verbose)

% Add somePath and its subfolders to the Matlab search path.
% A manual path save should be issued to save the path for future Matlab
% sessions.

pathStr = fvgenpath(somePath);

addpath(pathStr,'-begin');
sepIndices = find(pathStr==pathsep);

if verbose
   disp('****************************************************')
   disp('Added following directories to path:')
   startindex = 1;
   componentLength = diff([0 sepIndices]) - 1;
   for i=1:length(sepIndices) % -1 %commentend on 2010/04/01
       endIndex = startindex - 1 + componentLength(i);
       disp(pathStr(startindex:endIndex));
       startindex = endIndex + 2*length(pathsep);
   end
   disp(' ')
   disp('Save path manually for use in future Matlab sessions')
   disp('****************************************************')
end

end