function retp = center(s, fill)
    slen=length(s)+2;
    spad=' ';
    i=1;
    while i <= (80-slen)/2
       spad = [spad,fill];
       i=i+1;
    end
    sfmt=[spad,' ',s,' ',spad];
    retp = sfmt;
end % center %