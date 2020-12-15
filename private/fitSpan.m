function yOut = fitSpan(slope, domain)

yOut = slope*(domain(3:end) - domain(1)) + ...
    repmat(domain(2), [size(domain, 1)-2, 1]);