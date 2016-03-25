K = 10;
alpha_0 = 1.0;

n_max = 1000;

gamma   = 1:K;
delta = zeros(1,K-1);
alpha = alpha_0 .* ones(1,K);

%sample abundance (perf observations)
n_true = poissrnd(gamma);
y      = n_true;

%compute the true likelihood
ll_base = prod(poisspdf(n_true, gamma));

[~,~,ll_fa] = forward_messages(y, gamma, delta, alpha, n_max);
ll_fa = exp(ll_fa);

ll_gffa = gf_forward(y, gamma, alpha_0, delta);