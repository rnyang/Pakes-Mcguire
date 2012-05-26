function retp = getint(prompt, nlow, nhigh, default)
   load pmg.mat
   n=NOCHECK-1;
   while 1
%      print prompt$+' ';; %
      n=getn(prompt,default);
      if n~= floor(n)
         sprintf('Please re-enter, this parameter must be integer-valued');
         continue;
      end
      if (n<nlow && nlow~=NOCHECK) || (n>nhigh && nhigh~=NOCHECK)
         sprintf('Please re-enter - valid range for this parameter is ');
         if (nhigh==NOCHECK)
            sprintf('>');
            compact(nlow);
         elseif (nlow==NOCHECK)
            sprintf('<');
            compact(nhigh);
         else
            sprintf('%f - ',compact(nlow));
            compact(nhigh);
         end
      else
         retp = n;
         return % added to prevent infinite loop?
      end
   end
end % getint %