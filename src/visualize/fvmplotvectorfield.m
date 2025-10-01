function gh = fvmplotvectorfield(a,lw,c)
%FVMPLOTVECTORFIELD plots a vector field
%  
%  FVMPLOTVECTORFIELD(MYFIELD,LW,C)
%     MYFIELD is a vector field
%     LW is the line width in pts of the vectors
%     C specifies which component to plot.
%


   if ~isa(a,'Field')
      error('A Field must be supplied to create a vectorfield plot.')
   end

   if ~isvector(a)
      error('A Vector Field must be supplied to create a vectorfield plot.')
   end
   
   if ~haselements(a)
      error('The supplied field is empty. Cannot make the vectorfield plot.')
   end
   
   
   tmp = get(a);
   if nargin > 2 && ~isempty(c)
      vecdata = zeros(size(tmp));
      vecdata(c,:) = tmp(c,:);
      clear tmp
   else
      vecdata = tmp;
   end
   
%    zone = getzone(a);      
%    coords = getcoords(zone);
   coords = getcoords(a);
   
   gh = quiver(gca,coords(1,:),coords(2,:),vecdata(1,:),vecdata(2,:));   
   
   set(gh,'Color','k',...
          'LineWidth',lw,...
          'MaxHeadSize',0.08,...
          'AutoScaleFactor',0.75);
   

end