K = 10;
gamma_0 = 8;
alpha_0 = .5;

n_max = 1000;

gamma = [gamma_0, zeros(1,K-1)];
delta = ones(1,K-1);
alpha = alpha_0 .* ones(1,K);

%sample observations
n_true = poissrnd(gamma_0);
y = binornd(n_true, alpha_0, 1, K);

%compute the true likelihood
ll = zeros(n_max+1,1);
for n = 0:n_max
	ll(n+1) = poisspdf(n, gamma_0) * prod(binopdf(y, n, alpha_0));
end
ll_base = sum(ll);

[~,~,ll_fa] = forward_messages(y, gamma, delta, alpha, n_max);
ll_fa = exp(ll_fa);

ll_gffa = gf_forward(y, gamma, alpha_0, delta);