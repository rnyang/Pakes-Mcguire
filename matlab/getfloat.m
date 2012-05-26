function retp = getfloat(prompt, nlow, nhigh, default)
   load pmg.mat
   n=NOCHECK-1;
   while 1 % CREATE AN INFINITE LOOP?!
%      [prompt,' '] %
      n=getn(prompt, default);
      if (n<nlow && nlow~=NOCHECK) || (n>nhigh && nhigh~=NOCHECK);
         sprintf('Please re-enter - valid range for this parameter is ');
         if (nhigh==NOCHECK)
            sprintf('>');
            compact(nlow);
         elseif (nlow==NOCHECK)
            sprintf('<');
            compact(nhigh);
         else
            sprintf(' %f - ',compact(nlow));
            compact(nhigh);
         end
      else
         retp = n;
         return % maybe add a return here to prevent infinite loop?
      end
   end
end % getfloat %