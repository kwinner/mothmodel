function [ n ] = abundance( Q )
%compute how many individuals are alive between each interval

K = size(Q,1) - 1;

%each element of n is a 2D sum of a submatrix of Q
n = arrayfun(@(k) sum(sum(Q(1:k, k+1:K+1))), 1:K);

end

