function [y, n, S, Z, Q, eq_class] = samplestate( mu, sigma, lambda, N, alpha, T, varargin )

DEFAULT_OBS_DIST = 'bino';
DEFAULT_NSAMPLES = 1;

parser = inputParser;
addParamValue(parser, 'obs_dist', DEFAULT_OBS_DIST);
addParamValue(parser, 'nSamples', DEFAULT_NSAMPLES);

parser.parse(varargin{:});
obs_dist = parser.Results.obs_dist;
nSamples = parser.Results.nSamples;

%K = number of observations
K = numel(T);

%sample N birthtimes from Normal(mu, sigma)
S = normrnd(mu, sigma, N, nSamples);

%sample N lifespans from Exp(lambda)
Z = exprnd(lambda, N, nSamples);

%compute death times
D = S + Z;

%compute abundance
n = arrayfun(@(t) sum(S <= t & D >= t), T, 'UniformOutput', false);
n = cat(1, n{:})';

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
	Q = zeros(K+1, K+1, nSamples);
	%compute Q
	for iSample = 1:nSamples
		Q(:,:,iSample) = histc_2d(S(:,iSample), S(:, iSample) + Z(:, iSample), [-inf T inf], [-inf T inf]);
	end
	Q = squeeze(Q);
elseif nargout == 6 %special case to skip comptuting eq_class (which are only needed for some displays)
	Q = zeros(K+1, K+1, nSamples);
	eq_class = zeros(N, nSamples);
	for iSample = 1:nSamples
		%compute Q and eq_class
		[Q(:,:,iSample), ind_x, ind_y] = histc_2d(S(:, iSample), S(:, iSample) + Z(:, iSample), [-inf T inf], [-inf T inf]);

		%convert the x,y indices into equivalence classes
		[~,~,eq_class(:,iSample)] = unique(sub2ind([K+1, K+1], ind_x, ind_y));
	end
	Q = squeeze(Q);
	eq_class = squeeze(eq_class);
end

end

