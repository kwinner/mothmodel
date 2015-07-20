function [S, Z] = sample_pop( arrivalDistn, serviceDistn, N )
% SAMPLE_POP := Generate N individuals from the arrival and service time distributions
% [S, Z] = sample_pop( arrivalDistn, serviceDistn )
%          (sample a single individual)
% [S, Z] = sample_pop( arrivalDistn, serviceDistn, N )
%          (sample N individuals)
%    arrivalDistn  = a distribution object for the birth process (typically normal)
%    serviceDistn  = a distribution object for the death process (typically exponential)
%    N             = the total number of individuals (how many samples to draw)
%
%    S             = vector [N x 1] of individual birth times
%    Z             = vector [N x 1] of individual death times

if ~exist('N', 'var')
	N = 1;
end

S = arrivalDistn.random(N, 1);
Z = serviceDistn.random(N, 1);

end

