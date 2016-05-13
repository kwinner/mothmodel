K = 5;
alpha_0 = 0.1;

n_max = 50;

gamma   = 1:K;
delta = 0.5*ones(1,K-1);
alpha = alpha_0 .* ones(1,K);

% Observations
y = 10*ones(1,K);

% Run forward algorithm
[ll_gffa, a1, b1, f1, msg] = gf_forward(y, gamma, alpha_0, delta);
fprintf('Forward algorithm log-likelihood: %.4f\n', ll_gffa);

% Get marginals by performing tail elimination starting with the forward 
% message for all i=1:K

marg(K) = struct('a', 0, 'b', 0, 'f', 1);
for i=1:K
    m = struct();
    [ll_elim, m.a, m.b, m.f] = gf_tail_eliminate(y, gamma, alpha_0, delta, i, msg(i).a, msg(i).b, msg(i).f);
    marg(i) = m;
    
    fprintf('Marginal %d log-likelihood: %.4f\n', i, ll_elim);
end
