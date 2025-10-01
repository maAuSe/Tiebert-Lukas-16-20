function fvmplotfacenumbers(dom,fs)

fCoord = dom.fCoord;

for jF = 1:dom.nF
      
   text(fCoord(1,jF),fCoord(2,jF),num2str(jF),...
        'HorizontalAlignment','center',...
        'HitTest','off',...
        'FontAngle','normal',...
        'FontUnits','points',...        
        'FontSize',fs);
   
end


end