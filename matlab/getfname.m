function retp = getfname(prompt)
   sprintf('%s ',prompt);
   s=input('','s');
   s=lower(s);
   s=s(1,8);
   n=strfind(s,'.');
   n = n(1);
   if (n~=0)
       s=s(1,n-1);
   end
   sprintf(' ');
   retp = s;
end % getfname %