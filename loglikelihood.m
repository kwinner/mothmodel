function [ LL ] = loglikelihood (state, y, params)

LL = gammaln(params.N+1) + sum(alogb(state.q(:), state.p(:))) - sum(gammaln(state.q(:)+1)) + sum(dbinom(y, state.n, params.alpha, true));

% if isinf(LL) && LL < 0
%     LL = -10^10;
% end

end

function s = alogb(a, b)
s = a.*log(b);
s(a==0) = 0;
end