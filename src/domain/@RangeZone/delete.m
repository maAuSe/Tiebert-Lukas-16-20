function delete(This)

   dom = getdomain(This);
   % dom might be cleared earlier than the zones by the ML engine
   if isvalid(dom)
      undefinezone(dom,getid(This));
   end

end