function FieldObj = restrictto(This,zoneObj)
   
   thisZone = This.zone;
   TypeChecker.requireclass(zoneObj,class(thisZone));
   
   range  = thisZone.range;
   srange = zoneObj.range;
   if srange(1) < range(1) || srange(2) > range(2)
      error('Supplied zone must be a subrange.')
   end
   
   dsr1 = srange(1) - range(1) + 1;
   dsr2 = srange(2) - range(1) + 1;

   data = This.data;
   
   FieldObj = Field(zoneObj,This.rank);
   if ~isempty(data)
      set(FieldObj,data(:,dsr1:dsr2));
   end


end