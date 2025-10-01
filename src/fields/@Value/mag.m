function c = mag(a)

   dim = a.dim;
   
   if a.isscalar %---------------------------------------------------------
      % mag(s)
      datac = mag_s(get(a),getelcount(a));
      c = Value(0,dim);
      set(c,datac);
      return
   end
   
   if a.isvector %---------------------------------------------------------
      % mag(v)
      datac = mag_v(get(a),getelcount(a),dim);
      c = Value(0,dim);
      set(c,datac);
      return
   end
   
   if a.istensor %---------------------------------------------------------
      % mag(t)
      datac = mag_t(get(a),getelcount(a),dim);
      c = Value(0,dim);
      set(c,datac);
      return
   end

end