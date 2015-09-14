function y = logmnpdf( x, P )

y = gammaln(sum(x)+1) + sum(alogb(x(:), P(:))) - sum(gammaln(x(:)+1));

end


function s = alogb(a, b)
s = a.*log(b);
s(a==0) = 0;
end

