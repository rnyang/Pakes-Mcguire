function retp = mse(darray)
  % retp(padl(mean(darray),1,2)$+' ('$+padl(std(darray),1,2)$+')'); %  % RLG added std.errors %
  retp = [padr(mean(darray),8,3) , '  std: ' , padr(std(darray),8,3) , ' (' , padr(stds(darray)/sqrt(rows(darray)),8,5) , ')  range: ' , padr(min(darray),8,3) , ' ' , padr(max(darray),8,3)];
end