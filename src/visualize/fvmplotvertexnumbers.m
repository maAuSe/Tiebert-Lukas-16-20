function fvmplotvertexnumbers(dom,fs)

vCoord = dom.vCoord;

for jV = 1:dom.nV
      
   text(vCoord(1,jV),vCoord(2,jV),num2str(jV),...
        'HorizontalAlignment','center',...
        'HitTest','off',...
        'FontAngle','italic',...
        'FontUnits','points',...        
        'FontSize',fs);
   
end


end