function [A,b] = to_msparse(this)

   A = sparse(this.ir,this.jc,this.adata,...
              cast(this.n,'double'),cast(this.n,'double'));
   b = this.bdata;
   
end