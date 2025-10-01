function time(tinterval,annotationstring)

	fprintf('\t%7.4f s : %s\n\n',...
      round(10000*tinterval)/10000,annotationstring);

end