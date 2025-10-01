function c = ctranspose(a)

   tmp = get(a);
   if a.istensor
      tmp = transp_t(get(a),getelcount(a),a.dim);
   end
   c = Field(getzone(a),a.rank);
   set(c,tmp);

end