function c = ctranspose(a)

   tmp = get(a);
   if a.istensor
      tmp = transp_t(get(a),getelcount(a),a.dim);
   end
   c = Value(a.rank,a.dim);   
   set(c,tmp);

end