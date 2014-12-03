function [ LL ] = loglikelihood (state, y, params)

LL = gammaln(params.N+1) + sum(alogb(state.q(:), state.p(:))) - sum(gammaln(state.q(:)+1)) + sum(logbinopdf(y, state.n, params.alpha));

end

function s = alogb(a, b)
	s = a.*log(b);
	s(a==0) = 0;
end