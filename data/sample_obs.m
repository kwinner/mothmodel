function y = sample_obs( alpha, n, varargin )
% SAMPLE_OBS := Sample the observed counts from the abundance at each sampling time
% y = sample_obs( alpha, n )
% y = sample_obs( alpha, n, distn )
%
% INPUTS
% required:
%    alpha    = the observation probability
%    n        = vector [1 x K] of abundance at each observation time
%    distn    = the name of the distribution to use (Binomial or Poisson)
%
% OUTPUTS
%    y        = vector [1 x K] of observed counts at each observation time

DEFAULT_DISTN  = 'poisson';
POSSIBLE_DISTN = {'poisson', 'binomial'};

parser = inputParser;
addRequired(parser, 'alpha', @(x) isnumeric(x) && (0 <= x <= 1));
addRequired(parser, 'n',     @(x) isnumeric(x));
addOptional(parser, 'distn', DEFAULT_DISTN, @isstr);

%check inputs
parser.parse(alpha, n, varargin{:});
distn = parser.Results.distn;

switch validatestring(distn, POSSIBLE_DISTN, mfilename, 'distn')
case 'binomial'
	y = binornd(n, alpha);
case 'poisson'
	y = poissrnd(alpha .* n);
end

end

