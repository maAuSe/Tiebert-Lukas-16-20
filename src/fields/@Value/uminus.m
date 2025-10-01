function c = uminus(a)

   c = Value(a.rank,a.dim);
   set(c,-get(a));

end