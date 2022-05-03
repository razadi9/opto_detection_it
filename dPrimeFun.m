function d = dPrimeFun(IHit,IFA)
e= true;
d = nan;
if isempty(IHit )
    warning ('IHit is empty')
    e = false;
end
if isempty(IFA )
    warning ('IFA is empty')
    e = false;
end
if e
rateHit = nnz(IHit)/numel(IHit);
rateFA = nnz(IFA)/numel(IFA);

if rateHit == 1
    rateHit = .99;
elseif rateHit == 0
    rateHit = .01;
end

if rateFA == 1
    rateFA = .99;
elseif rateFA == 0
    rateFA = .01;
end
d= norminv(rateHit) - norminv(rateFA);
end

end