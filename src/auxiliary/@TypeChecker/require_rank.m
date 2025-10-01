function require_rank(var,allowedranks)

   % PRE: isa(var,'PhysicalData')

   rank = var.rank;
   if isempty(find(allowedranks==rank,1))
      supprankstr = sprintf('%d ',allowedranks);
      error('FVMLAB:UnsupportedRank',...
            'The rank of "%s" must be one of the following: %s.\n',...
            inputname(1),supprankstr);
   end
        
end