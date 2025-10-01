function implementationnotstarted()

    if true
        ME = MException('FvmLab:ImplementationNotStarted',...
               ['\nFVMLab:ImplementationNotStarted',...
                '\nImplementation of the required routine is not started.'...
                '\n']);
        throwAsCaller(ME);
    end

end