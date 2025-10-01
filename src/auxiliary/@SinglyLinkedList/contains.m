function flag = contains(this,Data_)

    Finger = this.Head;
    while ~isempty(Finger) && ~eq(Finger.Data,Data_)
        % Note that we check standard for reference equality on classes and
        % value equality on built-in data types. The class of Data might
        % overload the eq method to check for data equality. Its the class'
        % responsability!
        Finger = Finger.Next;        
    end % loop until node is found
    flag = ~isempty(Finger);

end % SinglyLinkedList.contains