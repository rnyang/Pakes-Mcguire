function retp = getnum2(prompt, default) % NEED TO REVISIT THIS FUNCTION
% con(1,1) replacement for UNIX %
   sprintf('%s : ',prompt);
   s=input('','s');
   if length(s)==0
      sprintf('  >> %s', compact(default));
      retp = default;
      return;
   end
   for i = 1:length(s)
       if sum(s(i)==['0','1','2','3','4','5','6','7','8','9','.','-','+']) ~= 1
           disp '\nContains non-number, please re-enter');
           retp = getnum2(prompt,default);
           return;
       end
   end     
   retp = str2double(s);
end % getnum2 %