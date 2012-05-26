function retp = compact(n)
% convert numeric scalar to C-style compact form %
%    sg=ftos(n,'%-*.*lG',16,8); 
   sg = num2str(n,6);
   sn='';
   i=1;
   while i<=length(sg)
      if ~strcmp(sg(i,1),' ')
         sn=[sn,sg(i,1)];
      end
      i=i+1;
   end
   retp = sn;
end % compact %