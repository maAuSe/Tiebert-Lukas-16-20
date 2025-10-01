function launchfvframework4students()
   % Set the path for this session

   type('license.txt');
   fprintf('\n');
   
   subdirs_to_add = {'src'};

   lastdir = pwd; % Save current path for later

   % Extract install directory (== this directory)
   fvmlabroot = mfilename('fullpath');
   endIndex = find(fvmlabroot==filesep,1,'last') - 1;
   fvmlabroot = fvmlabroot(1:endIndex);


   % Go to the install directory
   cd([fvmlabroot filesep 'src' filesep 'tools' filesep 'install'])
   % Run the addtopath script there with required path
   for ii=1:length(subdirs_to_add)
      addtopath([fvmlabroot filesep subdirs_to_add{ii}],false);
   end

    
   cd(lastdir) % Return to current path
   
   cd work

   clear fvmlabroot endIndex lastdir subdirs_to_add ii

fprintf('Paths successfully set.\n\n')

end