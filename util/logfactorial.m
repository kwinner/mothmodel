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
        if x == 0
            y = 0;
        else
            y = log(x) + gammaln(x);
        end
end

end

