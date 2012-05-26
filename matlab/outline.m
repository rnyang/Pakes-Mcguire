function retp = outline(title)
s=['---[',title,']'];
s=[s,strfill('-',79-length(s))];
retp = s;
end % outline %