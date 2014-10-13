function [ n ] = abundancy( q )

T = size(q,1) - 1;

%each element of n is a 2D sum of a submatrix of q
n = arrayfun(@(k) sum(sum(q(1:k, k+1:T+1))), 1:T);

%eliminate rounding error so that n aligns with y when alpha = 1
%generally, this rounding error shouldn't be a problem, but it could accumulate
% n = round(n);

end

