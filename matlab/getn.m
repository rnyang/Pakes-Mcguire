function retp = getn(prompt, default)
% con(1,1) replacement %
   load pmg.mat
   if ostype==DOS
      retp = getnum(prompt, default);
   else
      retp = getnum2(prompt, default);
   end
end % getn %