function c = magsq(a)

   dim = a.dim;
   
   if a.isscalar %---------------------------------------------------------
      % magsq(s)
      datac = magsq_s(get(a),getelcount(a));
      c = Value(0,dim);
      set(c,datac);
      return
   end
   
   if a.isvector %---------------------------------------------------------
      % magsq(v)
      datac = magsq_v(get(a),getelcount(a),dim);
      c = Value(0,dim);
      set(c,datac);
      return
   end
   
   if a.istensor %---------------------------------------------------------
      % magsq(t)
      datac = magsq_t(get(a),getelcount(a),dim);
      c = Value(0,dim);
      set(c,datac);
      return
   end

end