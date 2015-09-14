function n = abundance( S, Z, T )
% ABUNDANCE := Compute the number of individuals alive at each observation
% n = abundance( S, Z, T )
%
% INPUTS
% required
%    S = vector [N x 1] of individual birth times
%    Z = vector [N x 1] of individual death times
%    T = vector [1 x K] of observation times (sample times)
%
% OUTPUTS
%    n = vector [1 x K] of abundance
%
% note: this function is overloaded in the LV model, computing n from Q

assert(size(S, 1) == size(Z, 1));

%readability definitions
D = S + Z; %D = individual death times

%each n is the sum over the indicator function below for all t_k in T
n = arrayfun(@(t_k) sum(S <= t_k & D >= t_k), ...
	         T, ...
	         'UniformOutput', false);
n = cat(1, n{:})'; %concatenate and reshape the output to [1 x K]

end

