mu = 8; 
sigma = 4; 
lambda = 3;

alpha = 0.5;
N_hat = 80;

T = 1:4:20;
K = numel(T);

arrivalDistn = makedist('Normal', 'mu', mu, 'sigma', sigma); 
rateFunc = @arrivalDistn.pdf;

serviceDistn = makedist('Exp', 'mu', lambda);

gamma = immigration_rate(rateFunc, serviceDistn, T, N_hat);
delta = survival_prob(serviceDistn, T);

n_true = poissrnd(gamma); 
y = binornd(n_true, alpha, 1, K);