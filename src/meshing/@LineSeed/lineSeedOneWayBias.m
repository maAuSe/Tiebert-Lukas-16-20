function obj = lineSeedOneWayBias(originV,displV,...
            nParts,expFac,biasStart)
         
   obj = LineSeed();
   obj.nVertices = int32(nParts + 1);
   obj.nSegm = int32(nParts);
   obj.dX = zeros(obj.nVertices,1);
   obj.dY = zeros(obj.nVertices,1);
   obj.originX = originV(1);
   obj.originY = originV(2);

   expFrac = zeros(nParts+1,1);

   switch biasStart
      case 'o' %origin

         if expFac == 1
            cumExpLength = nParts;
            for ii=0:nParts
               expFrac(ii+1) = ii/cumExpLength;
            end
         else
            cumExpLength = (1 - expFac^nParts)/(1 - expFac);
            for ii=0:nParts
               expFrac(ii+1) = ((1 - expFac^ii)/(1 - expFac))/cumExpLength;
            end
         end

      case 'd' % destination

         if expFac == 1
            cumExpLength = nParts;
            for ii=1:nParts
               expFrac(ii) = ii/cumExpLength;
            end
         else
            cumExpLength = (1 - expFac^nParts)/(1 - expFac);
            for ii=0:nParts
               expFrac(nParts+1-ii) = 1 - ((1 - expFac^ii)/(1 - expFac))/cumExpLength;
            end                  
         end               

      otherwise
         error(['Specify ''o'' or ''d'' for bias start']); %#ok<NBRAK>
   end

   for ii = 0:nParts
      obj.dX(ii+1) = obj.dX(1) + displV(1)*expFrac(ii+1);
      obj.dY(ii+1) = obj.dY(1) + displV(2)*expFrac(ii+1);
   end

   obj.displX = obj.dX(end) - obj.dX(1);
   obj.displY = obj.dY(end) - obj.dY(1);         

   % % OLD code: introduced small but annoying errors.
   % % 'origin'
   % obj.dX(1) = originV(1);
   % obj.dY(1) = originV(2);
   % for i = 2:1:nParts
   %    obj.dX(i) = obj.dX(i-1) + L1x*expFac^(i-2);
   %    obj.dY(i) = obj.dY(i-1) + L1y*expFac^(i-2);
   % end
   % obj.dX(end) = originV(1) + displV(1);
   % obj.dY(end) = originV(2) + displV(2);

   % % OLD code: introduced small but annoying errors.
   % % 'destination'
   % obj.dX(end) = originV(1) + displV(1);
   % obj.dY(end) = originV(2) + displV(2);
   % for i = nParts:-1:2
   %    obj.dX(i) = obj.dX(i+1) - ...
   %       L1x*expansionFac^(nParts-i);
   %    obj.dY(i) = obj.dY(i+1) - ...
   %       L1y*expansionFac^(nParts-i);
   % end
   % obj.dX(1) = originV(1);
   % obj.dY(1) = originV(2);

end % lineSeedOneWayBias