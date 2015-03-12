function [S, Z] = sample_pop( arrivalDistn, arrivalParams, serviceDistn, serviceParams, N )
% SAMPLE_POP := Generate N individuals from the arrival and service time distributions
% [S, Z] = sample_pop( arrivalDistn, arrivalParams, serviceDistn, serviceParams, [N] )
%    arrivalDistn  = a distribution object for the birth process (typically normal)
%    arrivalParams = cell vector of parameters for the arrival distn
%    serviceDistn  = a distribution object for the death process (typically exponential)
%    serviceParams = cell vector of parameters for the death distn
%    N             = the total number of individuals (how many samples to draw)
%
%    S             = vector [N x 1] of individual birth times
%    Z             = vector [N x 1] of individual death times
%
% for additional details about the distribution object format, see README.txt

if ~exist(N, 'var')
	N = 1;
end

S = arrivalDistn.rand(arrivalParams{:}, N);
Z = serviceDistn.rand(serviceParams{:}, N);

end

