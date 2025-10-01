function c = magsq(a)

   dim = a.dim;
   
   if a.isscalar %---------------------------------------------------------
      % magsq(s)
      datac = magsq_s(get(a),getelcount(a));
      c = Field(getzone(a),0);
      set(c,datac);
      return
   end
   
   if a.isvector %---------------------------------------------------------
      % magsq(v)
      datac = magsq_v(get(a),getelcount(a),dim);
      c = Field(getzone(a),0);
      set(c,datac);
      return
   end
   
   if a.istensor %---------------------------------------------------------
      % magsq(t)
      datac = magsq_t(get(a),getelcount(a),dim);
      c = Field(getzone(a),0);
      set(c,datac);
      return
   end

end