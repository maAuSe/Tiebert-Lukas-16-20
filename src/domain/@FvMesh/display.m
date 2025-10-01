function display(This)

    fprintf('\n  "%s" with characteristics:',class(This))    
    fprintf('\n\t    dim : %g',This.dim);
    fprintf('\n\t     nC : %g',This.nC);    
    fprintf('\n\t     nF : %g',This.nF);
    fprintf('\n\t     nV : %g',This.nV);    
    fprintf('\n\t    nIf : %g',This.nIf);
    fprintf('\n\t    nBf : %g',This.nBf);
    fprintf('\n\n');
   
end