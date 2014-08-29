function [ Q, N, LL ] = mcmc( p, q_0, n_0, y, alpha, varargin )
%MCMC perform a full run of MCMC from q_0

DEFAULT_ITERATIONS = 1000;
DEFAULT_MOVES      = {'pairwise','shuffle','mergesplit'};

parser = inputParser;
addParamValue(parser, 'iterations', DEFAULT_ITERATIONS, @isnumeric);
addParamValue(parser, 'moves',      DEFAULT_MOVES,      @iscell);

parse(parser, varargin{:});

nIter = parser.Results.iterations;
moves = parser.Results.moves;

%% start sampling!
Q  = cell(1,nIter);  Q{1}  = q_0;
N  = cell(1,nIter);  N{1}  = n_0;
LL = zeros(1,nIter); LL(1) = loglikelihood(p,q_0,n_0,y,alpha);
i = 2;
while i <= nIter
    %sample
    [qprime, nprime, ~, success] = gibbsSample(p,Q{i-1},N{i-1},y,alpha,'moves',moves);
    
    if success
        %save result
        Q{i}  = qprime;
        N{i}  = nprime;
        LL(i) = loglikelihood(p,qprime,nprime,y,alpha);
        i = i + 1;
    else
        %discard sample (it was the same anyways) and don't count iteration
        continue
    end
end

end

