function [ y ] = logfactorial( x, mode )

if nargin < 2
    mode = 'gamma';
end

switch mode 
    case 'sumlog'
        y = sum(log(1:x));
    case 'stirling'
        %note: this is only an approximation
        y = x * log(x) - x + 1;
    case 'gamma'
        %note: gamma method is accurate, but undefined for x = 0
        y = log(x) + gammaln(x);
        y(x == 0) = 0;
end

end

