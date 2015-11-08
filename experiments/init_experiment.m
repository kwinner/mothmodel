%set default settings for moth model experiments

T     = 0:4:20;  %observation times
N_hat = 40;    %mean superpopulation size
n_max = N_hat; %cap on abundance at each time

alpha = 0.5;   %detection probability
mu      = 8; %mean arrival
sigma   = 4; %var of arrival
lambda  = 3; %mean service

%arrival rate function
arrivalDistn = makedist('Normal', 'mu', mu, 'sigma', sigma);
rateFunc     = @arrivalDistn.pdf;

%service distribution
serviceDistn = makedist('Exp', 'mu', lambda);