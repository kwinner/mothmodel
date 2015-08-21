%set default settings for moth model experiments

mu      = 8; %mean arrival
sigma   = 4; %var of arrival
lambda  = 3; %mean service

%arrival rate function
rateFunc = makedist('Normal', 'mu', mu, 'sigma', sigma);
rateFunc = @rateFunc.pdf;

%service distribution
serviceDistn = makedist('Exp', 'mu', lambda);

T     = 0:3:30;  %observation times
N_hat = 80;    %mean superpopulation size
n_max = N_hat; %cap on abundance at each time

alpha = 0.8;   %detection probability