function [ qprime, nprime, normalizedLikelihood ] = gibbs( p, q, n, y, alpha, a,b,c,d, varargin )

DEFAULT_LIKELIHOOD_MODE = 'simplified'; %simplified, full
DEFAULT_LOGSPACE        = false;

parser = inputParser;
addParamValue(parser,'likelihoodMode', DEFAULT_LIKELIHOOD_MODE, @(x) any(validatestring(x,{'simplified','full'})));
addParamValue(parser,'logspace',       DEFAULT_LOGSPACE,        @islogical);

parse(parser,varargin{:});

likelihoodMode = parser.Results.likelihoodMode;
logspace       = parser.Results.logspace;

T  = length(y);      %# of observations
NQ = (T+1)*(T+2)/2;  %size of upper triangular part of Q
N  = sum(sum(q));    %total # of individuals

%select two q to exchange
i = randsample(NQ,2,false);
i1 = i(1);
i2 = i(2);

%NOTE:
%  in all of the following code, we use the sense that individuals are
%  being transfered from outcome i1 to outcome i2, but since delta, the
%  number of individuals to transfer, can be negative, it is possible
%  that, for instance, the "positively" changed abundancies below are
%  simply those changed in the same sense as delta

%convert i1 and i2 from linear indices to subscripts
%[a,b] = ind2sub_triu(T+1, i1);
%[c,d] = ind2sub_triu(T+1, i2);

%find indices of abundancies which will have a "positive" change
nPos = union(a : (min(b,c) - 1), max(a,d) : (b - 1));

%find indices of abundancies which will have a "negative" change
nNeg = union(c : (min(d,a) - 1), max(c,b) : (d - 1));

%the y values place a minimum value on their corresponding n values
%the slack is how many individuals could be taken before y becomes impossible
%since the individuals are taken from every n, the tightest bound on any n is used
if isempty(nPos)
    nPosSlack = Inf;
else
    nPosSlack = min(n(nPos) - y(nPos));
end

if isempty(nNeg)
    nNegSlack = Inf;
else
    nNegSlack = min(n(nNeg) - y(nNeg));
end

deltaRange = max(-nPosSlack,-q(a,b)):min(nNegSlack,q(c,d));

if numel(deltaRange) <= 1
    qprime = q;
    nprime = n;
    normalizedLikelihood = mothLikelihood(p,q,n,y,alpha);
    return
end

%compute the likelihood for all possible deltas
%see notes for the reasoning behind each piece of this
if logspace
    unnormalizedLikelihood = zeros(1,numel(deltaRange));
else
    unnormalizedLikelihood =  ones(1,numel(deltaRange));
end
for iDelta = 1:numel(deltaRange)
    delta = deltaRange(iDelta);
    
    %will be 1 or 0 depending on UNLL init above
    prob = unnormalizedLikelihood(iDelta);
    
    switch likelihoodMode
        case 'simplified'
            if logspace
                %--simplified, logspace--
                prob = prob + delta * log(p(a,b)) - logfactorial(q(a,b)+delta);
                prob = prob - delta * log(p(c,d)) - logfactorial(q(c,d)-delta);
                for iPos = 1:numel(nPos)
                    n_i = n(nPos(iPos));
                    y_i = y(nPos(iPos));
                    
                    prob = prob + logfactorial(n_i + delta) - logfactorial(n_i + delta - y_i) + delta * log(1-alpha);
                end
                for iNeg = 1:numel(nNeg)
                    n_i = n(nNeg(iNeg));
                    y_i = y(nNeg(iNeg));
                    
                    prob = prob + logfactorial(n_i - delta) - logfactorial(n_i - delta - y_i) - delta * log(1-alpha);
                end
            else
                %--simplified, non-logspace--
                prob = prob * (p(a,b)^delta) / factorial(q(a,b)+delta);
                prob = prob * (p(c,d)^(-delta)) / factorial(q(c,d)-delta);
                prob = prob * (1-alpha)^(numel(nPos)*delta);
                prob = prob * (1-alpha)^(-numel(nNeg)*delta);
                for iPos = 1:numel(nPos)
                    n_i = n(nPos(iPos));
                    y_i = y(nPos(iPos));
                    
                    prob = prob * factorial(n_i + delta)/factorial(n_i + delta - y_i);
                end
                for iNeg = 1:numel(nNeg)
                    n_i = n(nNeg(iNeg));
                    y_i = y(nNeg(iNeg));
                    
                    prob = prob * factorial(n_i - delta)/factorial(n_i - delta - y_i);
                end
            end
        case 'full'
            %easiest to just calculate the modified q and n vectors
            qDelta = q;
            qDelta(a,b) = qDelta(a,b) + delta;
            qDelta(c,d) = qDelta(c,d) - delta;
            
            nDelta = abundancy(qDelta);
            
            if logspace
                %--full, logspace--
                prob = prob + logfactorial(N);
                for i = 1:numel(qDelta)
                    prob = prob + log(p(i)^qDelta(i)) - logfactorial(qDelta(i));
                end
                for k = 1:T
                    prob = prob + log(nchoosek(nDelta(k),y(k))) + y(k) * log(alpha) + (nDelta(k) - y(k)) * log(1-alpha);
                end
            else
                %--full, non-logspace--
                prob = prob * factorial(N);
                %linear indexing through all elements of q
                for i = 1:numel(qDelta)
                    prob = prob * (p(i) ^ qDelta(i)) / factorial(qDelta(i));
                end
                for k = 1:T
                    prob = prob * nchoosek(nDelta(k),y(k)) * alpha ^ y(k) * (1-alpha)^(nDelta(k)-y(k));
                end
            end            
    end
    
    unnormalizedLikelihood(iDelta) = prob;
end

%normalize the likelihood
if logspace
    shiftedUnnormalizedLikelihood = unnormalizedLikelihood - max(unnormalizedLikelihood);
    normalizedLikelihood = exp(shiftedUnnormalizedLikelihood) / sum(exp(shiftedUnnormalizedLikelihood));
else
    normalizedLikelihood = unnormalizedLikelihood / sum(unnormalizedLikelihood);
end

% if logspace && strcmp(likelihoodMode, 'full')
%     close all
%     plot(deltaRange, log(normalizedLikelihood))
%     pause
% end

% disp hello
%draw a delta according to the normalized likelihood
delta = deltaRange(randsample(numel(deltaRange),1,true,normalizedLikelihood));

%update q
qprime = q;
qprime(a,b) = qprime(a,b) + delta;
qprime(c,d) = qprime(c,d) - delta;

%update n
nprime = n;
nprime(nPos) = nprime(nPos) + delta;
nprime(nNeg) = nprime(nNeg) - delta;

end

