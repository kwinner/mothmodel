function [S, Z] = sample_pop( arrivalDistn, arrivalParams, serviceDistn, serviceParams, N )
% SAMPLE_POP := Generate N individuals from the arrival and service time distributions
% [S, Z] = sample_pop( arrivalDistn, arrivalParams, serviceDistn, serviceParams, [N] )
%    arrivalDistn  = a distribution object for the birth process (typically normal)
%    arrivalParams = rvector of parameters for the arrival distn
%    serviceDistn  = a distribution object for the death process (typically exponential)
%    serviceParams = rvector of parameters for the death distn
%    N             = the total number of individuals (how many samples to draw)
%
%    S             = cvector [N x 1] of individual birth times
%    Z             = cvector [N x 1] of individual death times

if ~exist(N, 'var')
	N = 1;
end

S = arrivalDistn.sample(N, deal(arrivalParams));
Z = serviceDistn.sample(N, deal(serviceParams));

end

