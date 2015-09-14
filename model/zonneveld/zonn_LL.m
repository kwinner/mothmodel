function LL = zonn_LL( y, arrivalDistn, serviceDistn, N_hat, alpha, T, varargin )
% ZONN_LL := compute the log likelihood of y given the Zonn model
% LL = zonn_LL( y, arrivalDistn, serviceDistn, T )
%
% INPUTS
% required:
%    y            = vector [1 x K] of observations of abundance
%    arrivalDistn = a distribution object for the birth process (typically normal)
%    serviceDistn = a distribution object for the death process (typically exponential)
%    N_hat        = mean (or actual) superpopulation size (positive int)
%    alpha        = detection probability (probability)
%    T            = vector [1 x K] of observation times (sample times)
%
% OUTPUTS
%    LL           = log likelihood of y

DEFAULT_DISTN  = 'poisson';
POSSIBLE_DISTN = {'poisson', 'binomial'};

parser = inputParser;
addOptional(parser, 'distn', DEFAULT_DISTN, @isstr)

parse(parser, varargin{:})
distn = parser.Results.distn;

%compute p_t, the prob that an indiv is alive at the times in T
p_t = presence_prob(arrivalDistn, serviceDistn, T);

%note, a binomial observation likelihood can be used, but is 0 if mean abundance is below a particular observation
%also, the binomial distn is undefined using the mean abundance and requires rounding...
switch validatestring(distn, POSSIBLE_DISTN, mfilename, 'distn')
case 'poisson'
	LL = sum(logpoisspdf(y, N_hat .* p_t .* alpha));
case 'binomial'
	LL = sum(logbinopdf(y, round(p_t .* N_hat), alpha));
end

end

