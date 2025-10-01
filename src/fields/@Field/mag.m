function c = mag(a)

   dim = a.dim;
   
   if a.isscalar %---------------------------------------------------------
      % mag(s)
      datac = mag_s(get(a),getelcount(a));
      c = Field(getzone(a),0);
      set(c,datac);
      return
   end
   
   if a.isvector %---------------------------------------------------------
      % mag(v)
      datac = mag_v(get(a),getelcount(a),dim);
      c = Field(getzone(a),0);
      set(c,datac);
      return
   end
   
   if a.istensor %---------------------------------------------------------
      % mag(t)
      datac = mag_t(get(a),getelcount(a),dim);
      c = Field(getzone(a),0);
      set(c,datac);
      return
   end

end