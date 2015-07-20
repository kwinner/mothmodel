function y = sample_obs( obsDistn, alpha, n )
% SAMPLE_OBS := Sample the observed counts from the abundance at each sampling time
% y = sample_obs( obsDistn, obsParams, n )
%    obsDistn = the name of the distribution to use (Binomial or Poisson)
%    alpha    = the observation probability
%    n        = vector [1 x K] of abundance at each observation time
%
%    y        = vector [1 x K] of observed counts at each observation time

parser = inputParser;
addRequired(parser, 'obsDistn', @(x) ismember(lower(x), {'binomial', 'poisson'}));
addRequired(parser, 'alpha',    @(x) isnumeric(x) && (0 <= x <= 1));
addRequired(parser, 'n',        @(x) isnumeric(x));

%check inputs
parser.parse(obsDistn, alpha, n);

switch(lower(obsDistn))
case 'binomial'
	y = binornd(n, alpha);
case 'poisson'
	y = poissrnd(alpha .* n);
end

end

