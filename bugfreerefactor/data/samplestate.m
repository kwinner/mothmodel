function [y, n, S, Z, Q, eq_class] = samplestate( mu, sigma, lambda, N, alpha, T, varargin )

DEFAULT_OBS_DIST = 'bino';

parser = inputParser;
addParamValue(parser, 'obs_dist', DEFAULT_OBS_DIST);

parser.parse(varargin{:});
obs_dist = parser.Results.obs_dist;

%K = number of observations
K = numel(T);

%sample N birthtimes from Normal(mu, sigma)
S = normrnd(mu, sigma, N, 1);

%sample N lifespans from Exp(lambda)
Z = exprnd(lambda, N, 1);

%compute death times
D = S + Z;

%compute abundance
n = arrayfun(@(t) sum(S <= t & D >= t), T);

%sample the observations
switch obs_dist
case 'bino'
	y = binornd(n, alpha);
case 'poiss'
	y = poissrnd(alpha .* n);
otherwise
	y = n;
end

if nargout == 5
	%compute Q
	Q = histc_2d(S, S+Z, [-inf T inf], [-inf T inf]);
elseif nargout == 6 %special case to skip comptuting eq_class (which are only needed for some displays)
	%compute Q and eq_class
	[Q, ind_x, ind_y] = histc_2d(S, S+Z, [-inf T inf], [-inf T inf]);

	%convert the x,y indices into equivalence classes
	[~,~,eq_class] = unique(sub2ind([K+1, K+1], ind_x, ind_y));
end

end

