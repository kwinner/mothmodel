function [ EQ, Q_history, runtime_history, success_history, attempts_history ] = mcmc( y, Q_0, P, N, alpha, T, varargin )

ALL_MOVES           = {'pair', 'shuffle', 'cycle', 'mergesplit'};

DEFAULT_NITERATIONS = 1000;      %# of mcmc iterations
DEFAULT_MOVES       = ALL_MOVES; %alternatively, take a subset of ALL_MOVES
DEFAULT_USE_ARS     = false;      %what LL sampler to use
DEFAULT_BURNIN      = 1;         %how many iterations to discard

parser = inputParser;
addParamValue(parser, 'nIterations', DEFAULT_NITERATIONS, @isnumeric);
addParamValue(parser, 'moves',       DEFAULT_MOVES,       @iscell);
addParamValue(parser, 'useARS',      DEFAULT_USE_ARS,     @islogical);
addParamValue(parser, 'burnin',      DEFAULT_BURNIN,      @isnumeric);

parse(parser, varargin{:});
nIterations = parser.Results.nIterations;
moves       = parser.Results.moves;
useARS      = parser.Results.useARS;
burnin      = parser.Results.burnin;

%initialize
Q_history    = cell(nIterations,1);
Q_history{1} = Q_0;
runtime_history = zeros(nIterations,1);
%runtime_history(1) = 0; %first Q was free
success_history = true(nIterations,1);
%success_history(1) = true;
attempts_history = zeros(nIterations,1);

for iteration = 2:nIterations
	disp(iteration)
	%do one sample
	starttime = tic;
	[Q_history{iteration}, success_history(iteration), attempts_history(iteration)] = gibbsSample(y, Q_history{iteration - 1}, P, N, alpha, T, moves, useARS);
	runtime_history(iteration) = toc(starttime);
end
EQ = mean(cat(3,Q_history{min(burnin+1, nIterations):end}),3);

end %/mcmc

%return Q_prime, a sampled state
%and success = was a valid move selected?
function [ Q_prime, success, attempts ] = gibbsSample( y, Q, P, N, alpha, T, moves, useARS )

K = length(y);      %# of observations
NQ = (K+1)*(K+2)/2; %size of upper triangular part of Q (# of valid bins)

success = false;
attempts = 1;

%compute abundancy of the current state
n = abundance(Q);

%select a move (uniformly at random)

move = struct;
move.type = moves{randi(numel(moves))};

while true

%select the cells to sample over
switch move.type
case 'pair'
	%sample two cells at random w/o replacement (using linear indexing into triu(Q))
	coords = randsample(NQ, 2, false);
	[move.a, move.b] = ind2sub_triu(K+1, coords(1));
	[move.c, move.d] = ind2sub_triu(K+1, coords(2));
case 'shuffle'
	%sample two cells from the diagonal of Q w/o replacement
	coords = randsample(K+1, 2, false);
	[move.a, move.b] = deal(coords(1));
	[move.c, move.d] = deal(coords(2));
case 'cycle'
	%sample 4 coordinates, with replacement
	coords = sort(randsample(K+1, 4, true));
	%we need to guarantee we selected 4 points
	%this is effectively just rejection sampling to avoid a pair/point selection
	while coords(1) == coords(2) || coords(3) == coords(4)
		coords = sort(randsample(K+1, 4, true));
	end

	%now that we have selected a valid 4 cycle...
	move.i = coords(1); move.iprime = coords(2);
	move.j = coords(3); move.jprime = coords(4);
case 'mergesplit'
	%basically the same as cycle above, but the second and third coordinates must be the same
	%this lands the bottom left corner of the cycle on the (unobserved) diagonal
	%note: we don't need to allow repeats, since repeating any coords would lead to a pair/point
	%      this also means the sample is guaranteed to work!
	coords = sort(randsample(K+1, 3, false));

	move.i = coords(1); move.iprime = coords(2);
	move.j = coords(2); move.jprime = coords(3);
end

%compute the range of legal values for delta
switch move.type
case {'pair', 'shuffle'}
	%compute how many individuals can be taken from Q(a,b)
	abSlack = Q(move.a, move.b);

	%find the minimum difference between the abundance and observations (which will be changed with this move)
	if move.a ~= move.b
		abImpactSet = union(move.a : (min(move.b, move.c) - 1), max(move.a, move.d) : (move.b - 1));
		abSlack = min([abSlack, n(abImpactSet) - y(abImpactSet)]);
	else
		abImpactSet = [];
	end

	%repeat for the cd move
	cdSlack = Q(move.c, move.d);
	if move.c ~= move.d
		cdImpactSet = union(move.c : (min(move.d, move.a) - 1), max(move.c, move.b) : (move.d - 1));
		cdSlack = min([cdSlack, n(cdImpactSet) - y(cdImpactSet)]);
	else
		cdImpactSet = [];
	end

	move.deltaRange = [-abSlack:cdSlack];

	%create masks for updating Q and n by delta
	move.updatemask_Q = sparse(K+1, K+1);
	move.updatemask_Q(move.a, move.b) = 1;
	move.updatemask_Q(move.c, move.d) = -1;

	move.updatemask_n = sparse(1,K);
	move.updatemask_n(abImpactSet) = 1;
	move.updatemask_n(cdImpactSet) = -1;
case {'cycle', 'mergesplit'}
	%slack in the top left/bottom right corners
	tlbrSlack = min(Q(move.i, move.j),      Q(move.iprime, move.jprime));
	%slack in the bottom left/top right corners
	bltrSlack = min(Q(move.i, move.jprime), Q(move.iprime, move.j     ));

	move.deltaRange = [-tlbrSlack:bltrSlack];

	%create masks for updating Q and n by delta
	move.updatemask_Q = sparse(K+1, K+1);
	move.updatemask_Q(move.i,      move.j)      = 1;
	move.updatemask_Q(move.iprime, move.jprime) = 1;
	move.updatemask_Q(move.i,      move.jprime) = -1;
	move.updatemask_Q(move.iprime, move.j)      = -1;

	%cycle moves don't change abundance...
	move.updatemask_n = sparse(1,K);
end

%no slack = no move to make, return same state
if numel(move.deltaRange) == 1
	Q_prime = Q;
	success = false;
	attempts = attempts + 1;
	continue
	% return
end

%construct our log likelihood function
logp = @(delta) latentvar_LL(Q + move.updatemask_Q .* delta, P) + latentvar_obs_LL(y, n + move.updatemask_n .* delta, alpha);

if useARS
	bounds = [min(move.deltaRange), max(move.deltaRange)];
	points = [];
	nsamples = 1;

	sampled_delta = discrete_ars(logp, bounds, points, nsamples);
else
	LL = arrayfun(@(delta) logp(delta), move.deltaRange);

	%this shouldn't ever happen if our methods are appropriately stable...
	% if any(LL == inf)
	% 	success = false;
	% end
	move.deltaRange = move.deltaRange(LL ~= inf);
	LL = LL(LL ~= inf);

	LL = normalizeLikelihood(LL);

	%LL should now be a probability vector, or stability failed somewhere...
	if any(LL < 0) || all(LL == 0) || any(isnan(LL))
		success = false;
		Q_prime = Q;
		attempts = attempts + 1;
		continue
		% return
	end

	sampled_delta = move.deltaRange(randsample(numel(move.deltaRange), 1, true, LL));
end

%apply the sampled_delta to get a new state
Q_prime = Q + sampled_delta .* move.updatemask_Q;

success = true;
break
end

end %/gibbsSample

function [ LL_prime ] = normalizeLikelihood( LL )

LL_prime = LL - max(LL); %shift for stability
LL_prime = exp(LL_prime);
LL_prime = LL_prime / sum(LL_prime);

end

