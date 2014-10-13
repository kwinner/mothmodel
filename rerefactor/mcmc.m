function [ state ] = mcmc( y, state_0, params, varargin )

DEFAULT_ITERATIONS = 1000;
DEFAULT_MOVES      = {'pairwise', 'shuffle', 'mergesplit'};

parser = inputParser;
addParamValue(parser, 'iterations', DEFAULT_ITERATIONS, @isnumeric);
addParamValue(parser, 'moves',      DEFAULT_MOVES,      @iscell);

parse(parser, varargin{:});
nIterations = parser.Results.iterations;
moves       = parser.Results.moves;

%initialize
state = state_0;

iter = 2;
while iter <= nIterations + 1
	%do one sample
	[state(iter), success] = gibbsSample(y, state(iter-1), params, 'moves', moves);

	%if the sample was successful, move on
	%unsuccessful moves are, for instance, moves which select an illegal pair to xfer
	%or which have no available slack
	if ~success
		continue
	else
		iter = iter + 1;
	end
end

end

function [ state, success ] = gibbsSample( y, state_0, params, varargin )

DEFAULT_MOVES = {'pairwise', 'shuffle', 'mergesplit'};

parser = inputParser;
addParamValue(parser, 'moves', DEFAULT_MOVES, @iscell);

parse(parser, varargin{:});
moves = parser.Results.moves;

T  = length(y);     %# of observations
NQ = (T+1)*(T+2)/2; %size of upper triangular part of Q (# of valid outcomes)

%select a move (uniformly at random)
move = moves{randi(numel(moves))};

switch move
	case 'pairwise'
		%select two cells in q to rebalance
		i = randsample(NQ, 2, false);
		move = struct;
		[move.a, move.b] = ind2sub_triu(T+1, i(1));
		[move.c, move.d] = ind2sub_triu(T+1, i(2));

		%%sample delta, the # of indivs to move from q(a,b) to q(c,d)

		%the slack computations below rely on the fact that if you xfer some
		%mass from q(a,b) to q(c,d), some observations increase, some decrease, and some see both changes
		%nAB/nCD are the cells that change an amount equal to the change in q(a,b)/q(c,d) respectively
		%and the corresponding slack is how far those abundancies are over the corresponding observations

		%how much slack is there in the abundancies that change with ab
		move.nAB = union(move.a : (min(move.b,move.c) - 1), max(move.a,move.d) : (move.b - 1));
		if isempty(move.nAB) nABSlack = Inf; %the cells that change with AB are a subset of those that change with CD
		else
			nABSlack = min(state_0.n(move.nAB) - y(move.nAB));
		end

		%how much slack is there in the abundancies that change with cd
		move.nCD = union(move.c : (min(move.d,move.a) - 1), max(move.c,move.b) : (move.d - 1));
		if isempty(move.nCD) nCDSlack = Inf; %the cells that change with CD are a subset of those that change with AB
		else
			nCDSlack = min(state_0.n(move.nCD) - y(move.nCD));
		end

		abMax = min(state_0.q(move.a,move.b), nABSlack);
		cdMax = min(state_0.q(move.c,move.d), nCDSlack);

		%if there is no slack available, then no move can be made
		if abMax == 0 && cdMax == 0
			success = false; state = state_0;
			return
		end

		%indices of P/Q that change with delta
		move.changedPQ = sub2ind([T+1,T+1], [move.a move.c], [move.b move.d]);
		move.changedNY = union(move.nAB, move.nCD);

		%range of values for delta
		%max negative is the maximum # that can be taken from ab
		%max positive is the maximum # that can be taken from cd
		deltaRange = [-abMax:cdMax];

		%compute the LL over all possible delta
		move.applydelta = @applydelta_pairwise;
		UNLL = arrayfun(@(delta) loglikelihoodofdelta(delta, state_0, y, params, move), deltaRange);

		%anywhere that UNLL = -inf is impossible and should be removed
		%this typically only occurs when params.alpha = 1 and delta makes the observations impossible
		deltaRange = deltaRange(UNLL ~= inf);
		UNLL = UNLL(UNLL ~= inf);

		NLL = normalizeLikelihood(UNLL);

		%randomly sample delta from the normalized LL
		delta = deltaRange(randsample(numel(deltaRange), 1, true, NLL));

		state = applydelta_pairwise(delta, state_0, move);
		success = true;
	case 'shuffle'
		%note: this move is very similar to 'pairwise' above, just with *many* simplifications
		%select two unobserved q to shuffle
		i = randsample(T, 2, false);
		move = struct;
		move.a = i(1); move.b = i(1);
		move.c = i(2); move.d = i(2);

		%since abundancy won't change w/ a shuffle, the delta limits are easy
		deltaRange = -state_0.q(move.a, move.b):state_0.q(move.c, move.d);

		%compute the LL over all possible delta
		move.applydelta = @applydelta_pairwise;
		UNLL = arrayfun(@(delta) loglikelihoodofdelta(delta, state_0, y, params, move), deltaRange);

		%unlike the pairwise case, this should never be triggered, but better safe than sorry!
		deltaRange = deltaRange(UNLL ~= inf);
		UNLL = UNLL(UNLL ~= inf);

		NLL = normalizeLikelihood(UNLL);

		%randomly sample delta from the normalized LL
		delta = deltaRange(randsample(numel(deltaRange), 1, true, NLL));

		state = applydelta_pairwise(delta, state_0, move);
		success = true;
	case 'mergesplit'
		%note: this is also an abundancy-preserving move, but with more complex
		%bound on delta
		%since it is slightly differeny from above,
		%positive delta in a m/s corresponds to increasing q(a,c)

		%sample three indices, i is the smallest, k is the largest
		move = struct;
		temp = sort(randsample(T+1, 3, false));
		move.i = temp(1); move.j = temp(2); move.k = temp(3);

		%smaller of the two bins being merged
		maxMerge = min(state_0.q(move.i,move.j), state_0.q(move.j,move.k));
		%smaller of the large bin, or the unobserved slack at j,j
		maxSplit = min(state_0.q(move.i,move.k), state_0.q(move.j,move.j));

		deltaRange = -maxSplit:maxMerge;

		%compute the LL over all possible delta
		move.applydelta = @applydelta_ms;
		UNLL = arrayfun(@(delta) loglikelihoodofdelta(delta, state_0, y, params, move), deltaRange);

		%this should also never happen w/ a ms
		deltaRange = deltaRange(UNLL ~= inf);
		UNLL = UNLL(UNLL ~= inf);

		NLL = normalizeLikelihood(UNLL);

		%randomly sample delta from the normalized LL
		delta = deltaRange(randsample(numel(deltaRange), 1, true, NLL));

		state = applydelta_ms(delta, state_0, move);
		success = true;
end

end

%%apply the move delta to state according to move
function [ state_prime ] = applydelta_pairwise( delta, state, move )

state_prime = state;
state_prime.q(move.a,move.b) = state.q(move.a,move.b) + delta;
state_prime.q(move.c,move.d) = state.q(move.c,move.d) - delta;

if isfield(move,'nAB')
	state_prime.n(move.nAB) = state.n(move.nAB) + delta;
end
if isfield(move,'nCD')
	state_prime.n(move.nCD) = state.n(move.nCD) - delta;
end

end

function [ state_prime ] = applydelta_ms( delta, state, move )

state_prime = state;
state_prime.q(move.i, move.k) = state.q(move.i, move.k) + delta;
state_prime.q(move.j, move.j) = state.q(move.j, move.j) + delta;
state_prime.q(move.i, move.j) = state.q(move.i, move.j) - delta;
state_prime.q(move.j, move.k) = state.q(move.j, move.k) - delta;

end

%%apply the move delta to state, then compute the UNLL
function [ UNLL ] = loglikelihoodofdelta( delta, state, y, params, move )

state_prime = move.applydelta(delta, state, move);

UNLL = loglikelihood(state_prime, y, params);

end

function [ NLL ] = normalizeLikelihood( UNLL )

NLL = UNLL - max(UNLL); %shift for stability
NLL = exp(NLL);
NLL = NLL / sum(NLL);

end