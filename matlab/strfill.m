function retp = strfill(fill, len)
    s=' ';
    i=1;
    while i <= len
       s = [s,fill];
       i=i+1;
    end
    retp = s;
end % strfill %