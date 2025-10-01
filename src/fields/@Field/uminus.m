function c = uminus(a)

   c = Field(getzone(a),a.rank);
   set(c,-get(a));

end