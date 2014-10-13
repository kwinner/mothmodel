function [ qprime, nprime, normalizedLikelihood, success, move, deltaRange ] = gibbsSample( p, q, n, y, alpha, varargin )
%GIBBSSAMPLE sample a new state qprime by taking one step from q
% this sampler can take three kinds of moves:
%  pairwise    - take any two outcomes, redistribute the individuals in those two
%  shuffle     - take two unobserved outcomes, redistribute the individuals
%  merge/split - take individuals from one long outcome and split them into individuals in 2 component outcomes, or vice versa

%% param processing
DEFAULT_MOVES = {'pairwise','shuffle','mergesplit'};

parser = inputParser;
addParamValue(parser,'moves',DEFAULT_MOVES,@iscell);

parse(parser,varargin{:});

MOVES = parser.Results.moves;

%% constants
T  = length(y);      %# of observations
NQ = (T+1)*(T+2)/2;  %size of upper triangular part of Q (# of valid outcomes)
N  = sum(sum(q));    %total # of individuals

%% select a move (uniformly at random)
move = MOVES{randi(numel(MOVES))};

%DEBUG
% move = 'mergesplit';

switch move
    case 'pairwise'
        %select two q to exchange
        i = randsample(NQ, 2, false);
        [a, b] = ind2sub_triu(T+1, i(1));
        [c, d] = ind2sub_triu(T+1, i(2));
        
        %compute the range of valid values for delta, the # of indivs to move from q(c,d) to q(a,b)
        posMax = fix(q(a,b));
        negMax = fix(q(c,d));
        
        %further constrain posMax and negMax by their impact on n
        %each n_i must remain higher than the corresponding y_i
        nPos = union(a : (min(b,c) - 1), max(a,d) : (b - 1));
        if isempty(nPos) nPosSlack = Inf;
        else
            nPosSlack = fix(min(n(nPos) - y(nPos))); %the smallest difference between an nPos and its yPos
        end
        
        nNeg = union(c : (min(d,a) - 1), max(c,b) : (d - 1));
        if isempty(nNeg) nNegSlack = Inf;
        else
            nNegSlack = fix(min(n(nNeg) - y(nNeg)));
        end
        
        %max and min = smallest (closest to 0) values for delta
        posMax = min(posMax,nPosSlack);
        negMax = min(negMax,nNegSlack);
        
        %posMax = the most that can be taken from q(a,b)
        %negMax = the most that can be taken from q(c,d)
        deltaRange = -posMax:negMax;
        
        %if there is only one (or zero...) possible delta, then delta is too constrained
        %and the move is a failure, so return without taking a step
        if numel(deltaRange) <= 1
            qprime = q;
            nprime = n;
            normalizedLikelihood = 1;
            success = false;
            return
        end
        
        %% otherwise, for each possible delta in deltaRange, compute the UNLL (unnormalized log likelihood)
        UNLL = zeros(size(deltaRange));
        qprime = q;
        nprime = n;
        changedPQ = sub2ind(size(p),[a c],[b d]);
        changedNY = union(nPos, nNeg);
        for iDelta = 1:numel(deltaRange)
            delta = deltaRange(iDelta);
            qprime(a,b) = q(a,b) + delta;
            qprime(c,d) = q(c,d) - delta;
            
            nprime(nPos) = n(nPos) + delta;
            nprime(nNeg) = n(nNeg) - delta;
            
            UNLL(iDelta) = loglikelihood(p(changedPQ), qprime(changedPQ), nprime(changedNY), y(changedNY), alpha);
        end

        %normalize UNLL
        normalizedLikelihood = normalizeLikelihood(UNLL);
        
        if ~any(normalizedLikelihood)
            true
        end
        
        %% pick a delta, and apply that, completing the move
        delta = deltaRange(randsample(numel(deltaRange),1,true,normalizedLikelihood));
        
        qprime = q;
        qprime(a,b) = q(a,b) + delta;
        qprime(c,d) = q(c,d) - delta;
        
        nprime = n;
        nprime(nPos) = n(nPos) + delta;
        nprime(nNeg) = n(nNeg) - delta;
        
        success = true;
        return
    case 'shuffle'
        %note: this case is very similar to 'pairwise' above, just with *many* simplifications
        %select two unobserved q to shuffle
        i = randsample(T, 2, false); 
        a = i(1);
        b = i(2);
        
        %since abundancy won't change w/ a shuffle, the delta limits are easy
        deltaRange = -q(a,a):q(b,b);
        
        %this can only happen if q(a,a) == q(b,b) == 0 in a shuffle
        if numel(deltaRange) <= 1
            qprime = q;
            nprime = n;
            normalizedLikelihood = 1;
            success = false;
            return
        end
        
        %% compute the UNLL for each possible delta
        UNLL = zeros(size(deltaRange));
        for iDelta = 1:numel(deltaRange)
            delta = deltaRange(iDelta);
            UNLL(iDelta) = loglikelihood([p(a,a) p(b,b)], [q(a,a) + delta, q(b,b) - delta], [], [], alpha);
        end
        
        %normalize the UNLL
        normalizedLikelihood = normalizeLikelihood(UNLL);
        
        %% pick a delta, and apply that, completing the move
        delta = deltaRange(randsample(numel(deltaRange),1,true,normalizedLikelihood));
        
        qprime = q;
        qprime(a,a) = qprime(a,a) + delta;
        qprime(b,b) = qprime(b,b) - delta;
        
        nprime = n; %no change to n
        
        success = true;
        return
    case 'mergesplit'
        %note: this is also an abundancy-preserving move, but with more complicated delta bounds
        %since it is slightly different from above, positive delta in a m/s here corresponds to increasing the individuals in the large outcome
        %sample three indices, i is the smallest, k is the largest
        temp = sort(randsample(T, 3, false));
        i = temp(1); j = temp(2); k = temp(3);
        
        maxMerge = min(q(i,j), q(j,k)); %smaller of the two bins that would be merged
        maxSplit = min(q(i,k), q(j,j)); %smaller of the large bin, or the unobserved indivs at j,j needed to make a split
        
        deltaRange = -maxSplit:maxMerge;
        
        %as usual, if there are no valid deltas, return a fail
        if numel(deltaRange) <= 1
            qprime = q;
            nprime = n;
            normalizedLikelihood = 1;
            success = false;
            return
        end
        
        %% for each delta, compute UNLL
        UNLL = zeros(size(deltaRange));
        for iDelta = 1:numel(deltaRange)
            delta = deltaRange(iDelta);
            UNLL(iDelta) = loglikelihood([p(i,j),         p(j,k),         p(i,k),         p(j,j)], ...
                                         [q(i,j) - delta, q(j,k) - delta, q(i,k) + delta, q(j,j) + delta], ...
                                         [], [], alpha);
        end
        
        %normalize the UNLL
        normalizedLikelihood = normalizeLikelihood(UNLL);
        
        if ~any(normalizedLikelihood)
            true
        end
        
        %% pick a delta, and apply that, completing the move
        delta = deltaRange(randsample(numel(deltaRange),1,true,normalizedLikelihood));
        
        qprime = q;
        qprime(i,j) = qprime(i,j) - delta;
        qprime(j,k) = qprime(j,k) - delta;
        qprime(i,k) = qprime(i,k) + delta;
        qprime(j,j) = qprime(j,j) + delta;
        
        nprime = n; %no change to n
        
        success = true;
end

end

function [ NLL ] = normalizeLikelihood( UNLL )

NLL = UNLL - max(UNLL); %shift for stability
NLL = exp(NLL);
NLL = NLL / sum(NLL);

end