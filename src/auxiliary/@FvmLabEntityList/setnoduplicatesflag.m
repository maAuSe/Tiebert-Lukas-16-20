function setnoduplicatesflag(this,flag)

    TypeChecker.requireclass(flag,'logical');
    TypeChecker.requiresingleelem(flag);
    
    this.noDuplicatesFlag = flag;

end