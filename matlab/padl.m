function retp = padl(n,w,p)
% n left-padded to w with precision p %
%    sg=ftos(n,'%*.*lf',w,p);
   sg = num2str(n);
%    sg = sg(1:w);
   sn='';
   i=1;
   while i<=length(sg)
      if ~strcmp(sg(i,1),' ');
         sn=[sn,sg(i,1)];
      end
      i=i+1;
   end
   while length(sn)<=w
      sn=[sn,' '];
   end
   retp = sn;
end % padl %