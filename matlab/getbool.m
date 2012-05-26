function retp = getbool(prompt)
   sprintf('%s (y/n) ? ',prompt);
   s=input('','s');
   retp = iif(strcmp(s,'y') || strcmp(s,'Y'), 1, 0);
end % getbool %