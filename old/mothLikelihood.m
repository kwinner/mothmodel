function [ LL ] = mothLikelihood( p, q, n, y, alpha )

%prob of q given p
%LL = log(mnpdf(reshape(q,1,numel(q)), reshape(p,1,numel(p))));
LL = mnpdf(reshape(q,1,numel(q)), reshape(p,1,numel(p)));

%prob of y given n, alpha
for i = 1:numel(y)
    LL = LL * binopdf(y(i), n(i), alpha);
end

end