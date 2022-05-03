function s = num2strComma(n)
if n ~= fix(n)
    error('only integer')
end
dec = 3;
c   = fix(log10(n)+1);
fmt = [repmat('%c',1,mod(c,dec)) repmat(',%c%c%c',1,fix(c/dec)) '%s'];
s = sprintf(fmt, sprintf('%d',n));
if s(1) == ','
    s(1) = [];
end
