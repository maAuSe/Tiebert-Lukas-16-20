function implementationnotfinished()

    if true
        ME = MException('FvmLab:ImplementationNotFinished',...
               ['\nFVMLab:ImplementationNotFinished',...
                '\nImplementation of the required routine is not finished.'...
                '\n']);
        throwAsCaller(ME);
    end

end