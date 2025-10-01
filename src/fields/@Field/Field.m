classdef Field < Value & IField
   
   
properties (GetAccess=public,SetAccess=protected)

   dom         % FvDomain
   % dom that holds the zone
   zone           % FvZone
   % zone on which this field is defined
   
   elcountzone    % int32
   % Holds the number of data elements in order to fill up the zone.
   
   % NOTE: we must overload the majority of the routines because the target
   % class changes, and hence the constructor changes. We might create
   % a sort of empty copy constructors with possibility to optionally alter
   % the rank and dimension. This way the base class is retained!!!!
   
end % PROPERTIES



methods (Access=public,Static=true)
   
   
   
end % STATIC METHODS



methods (Access=public)

   function Obj = Field(zone,rank)
   %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   % Field object constructor. Creates an empty container for storing
   % values of the given rank that represent a field on the given zone.
   % POST:
   %  A valid Value object iff DATA meets all the following conditions:   
   %  (1) rank is one of the following 0,1,2
   %  (2) zone is a FvZone
   %
   %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      Obj = Obj@Value(rank,getdim(zone));
      Obj = Obj@IField;
      
      Obj.dom  = zone.dom;
      Obj.zone = zone;
      Obj.elcountzone = getelcount(zone);

   end

   % Special Field methods
   
   % Implementation of methods in IField
   DomObj      = getdomain(This)
   ZoneObj     = getzone(This)
   coordsdata  = getcoords(This)
   n           = getelcountzone(This)
   FieldObj    = restrictto(This,zoneObj)
   
   % Overload of methods in Value
   % data = get(This,c1,c2)            % NOT overloaded
   % data = extract(This)              % NOT overloaded
   set(This,data,c1,c2)
   reset(This,eldata,c1,c2)
   c = comp(This,c1,c2)
   
   % c = plus(a,b)                     % NOT overloaded
   % c = minus(a,b)                    % NOT overloaded
   % c = mtimes(a,b)                   % NOT overloaded
   c = uminus(a)
   c = uplus(a)
   c = ctranspose(a)

   c = magsq(a)
   c = mag(a)
    
   % Implementation of auxiliary methods for binary operators
   c = lhs_plus(a,b)
   c = rhs_plus(b,a)
   c = lhs_minus(a,b)
   c = rhs_minus(b,a)
   c = lhs_mtimes(a,b)
   c = rhs_mtimes(b,a)
   c = lhs_dot(a,b)
   c = rhs_dot(b,a)
      
   % Implementation of methods PhysicalData
   % elcount = getelcount(This)        % NOT overloaded
   % haselflag = haselements(This)     % NOT overloaded


end % METHODS

   
   
   
end












