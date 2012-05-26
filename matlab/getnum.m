function retp = getnum(prompt, default)
% con(1,1) replacement for DOS %
   sprintf('%s : ',prompt);
   while ch ~= 13 % Enter %
      s=input('','s');
   end
   if length(s)==0
      sprintf('%s',compact(default));
      retp = default;
   else
      retp = str2double(s);
   end
end % getnum %