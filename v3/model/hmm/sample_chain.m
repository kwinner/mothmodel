function n = sample_chain( rateFunc, serviceDistn, T, N_hat )

n = zeros(size(T));

%sample the # of initial individuals
P_0 = initial_distn(rateFunc, serviceDistn, T(1), N_hat);
n(1) = randsample(0:(numel(P_0)-1),1,true,P_0);

%compute transitions
for k = 2:numel(T)
	tk1 = T(k-1);
	tk2 = T(k);

	%sample transition
	PT_k = transition_distn(rateFunc, serviceDistn, tk1, tk2, N_hat, n(k-1));
	n(k) = randsample(0:(numel(PT_k)-1),1,true,PT_k);
end

end

