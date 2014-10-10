function [ LL ] = loglikelihood( p, q, n, y, alpha, varargin )
%compute the LOG LIKELIHOOD
%note: p,q,n,y do not need to be complete. 
% you can use corresponding subsets to compute an unnormalized LL faster

%the calculation mode.
%  pdf is much quicker, but has numerical under/overflow problems in some cases
%  stable is more stable, but not optimized for speed
DEFAULT_MODE = 'stable'; %stable, pdf

parser = inputParser;
addParamValue(parser,'mode',DEFAULT_MODE, @(x) any(validatestring(x,{'stable','pdf'})));

parse(parser,varargin{:});

mode = parser.Results.mode;

switch mode
    case 'pdf'
        %prob of q given p
        LL = log(mnpdf(reshape(q,1,numel(q)), reshape(p,1,numel(p))));
        
        %prob of y given n, alpha
        for i = 1:numel(y)
            LL = LL + log(binopdf(y(i), n(i), alpha));
        end
    case 'stable'
        N = sum(q(:));
        
        LL = logfactorial(N);
        
        if ~isempty(p)
            %compute each prior
            %stable log of (p_i^q_i)/q_i!
            LL = LL + sum(arrayfun(@prior,p(p(:) ~= 0),q(p(:) ~= 0)));
        end
        
        if ~isempty(y)
            %compute each posterior
            %stable log of (n_i nchoosek y_i) * alpha^y_i * (1-alpha)^(n_i-y_i)
            LL = LL + sum(arrayfun(@(y_i,n_i) posterior(y_i,n_i,alpha),y(:),n(:)));
        end
end
end

%% compute the prior probability of p_i,q_i
function [ prob ] = prior( p_i, q_i )
if p_i == 0 && q_i == 0
    prob = 0;
else
    prob = q_i*log(p_i) - logfactorial(q_i);
end
end

%% compute the posterior probability of y_i,n_i
function [ prob ] = posterior( y_i, n_i, alpha )
epsilon = 1e-6;
%fix for rounding error
difference = n_i - y_i;
if abs(difference) < epsilon
    difference = 0;
end
prob = logfactorial(n_i) - logfactorial(y_i) - logfactorial(difference);
if ~(y_i == 0 && alpha == 0) && ~(y_i == n_i && alpha == 1)
    prob = prob + y_i * log(alpha) + (difference) * log(1-alpha);    
end
end