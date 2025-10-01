classdef StopfileMonitor
    
   
   properties
      
      stopfn = 'stop';
      dontstopfn = 'dontstop';
      
   end

   methods
      
      function This = StopfileMonitor()
         
         % clear existing files
         if exist(This.stopfn,'file')
            delete(This.stopfn);
         end
         if exist(This.dontstopfn,'file')
            delete(This.dontstopfn);
         end
         
         fid = fopen(This.dontstopfn,'w');
         fclose(fid);
         
      end
      
      
      function stoprequestflag = checkstoprequest(This)
         
         fid = fopen(This.stopfn,'r');
         if fid > 0
            fclose(fid);
            stoprequestflag = true;
         else
            stoprequestflag = false;
         end
         
      end
      
      
   end
    
        
end