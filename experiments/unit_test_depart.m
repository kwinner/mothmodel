K = 10;
gamma_0 = 8;
alpha_0 = 1.0;
delta_hat = 0.75;

n_max = 1000;

gamma = [gamma_0, zeros(1,K-1)];
delta = delta_hat .* ones(1,K-1);
alpha = alpha_0 .* ones(1,K);

%sample abundance (perf observations)
n_true(1) = poissrnd(gamma_0);
for i = 2:K
	n_true(i) = binornd(n_true(i-1), delta(i-1));
end
y = n_true;

%compute the true likelihood
ll_base = poisspdf(n_true(1), gamma_0) .* prod(binopdf(n_true(2:end), n_true(1:end-1), delta));

[~,~,ll_fa] = forward_messages(y, gamma, delta, alpha, n_max);
ll_fa = exp(ll_fa);

ll_gffa = gf_forward(y, gamma, alpha_0, delta);