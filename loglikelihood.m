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
            prior = @(p_i, q_i) q_i*log(p_i) - logfactorial(q_i);
            LL = LL + sum(arrayfun(prior,p(p(:) ~= 0),q(p(:) ~= 0)));
        end
        
        if ~isempty(y)
            %compute each posterior
            %stable log of (n_i nchoosek y_i) * alpha^y_i * (1-alpha)^(n_i-y_i)
            posterior = @(y_i, n_i) logfactorial(n_i) - logfactorial(y_i) - logfactorial(n_i-y_i) ...
                + y_i * log(alpha) + (n_i-y_i) * log(1-alpha);
            LL = LL + sum(arrayfun(posterior,y(:),n(:)));
        end
end
end

