function retp = iif(a,b,c) 
% analog of C (a) ? b : c statement %
   if a
      retp = b;
   else
      retp = c;
   end
end % iif %