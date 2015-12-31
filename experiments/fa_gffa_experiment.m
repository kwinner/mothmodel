T     = 0:4:20;  %observation times

alpha = 0.5;   %detection probability
mu      = 8; %mean arrival
sigma   = 4; %var of arrival
lambda  = 3; %mean service

N_hat_vals = 50:50:1000;
runtime_fa   = zeros(10,numel(N_hat_vals));
runtime_gffa = zeros(10,numel(N_hat_vals));
for iN_hat = 1:numel(N_hat_vals)
	N_hat = N_hat_vals(iN_hat);
	N_max = N_hat;
	%arrival rate function
	arrivalDistn = makedist('Normal', 'mu', mu, 'sigma', sigma);
	rateFunc     = @arrivalDistn.pdf;

	%service distribution
	serviceDistn = makedist('Exp', 'mu', lambda);

	N_hat = N_hat_vals(iN_hat);    %mean superpopulation size
	n_max = N_hat; %cap on abundance at each time

	for irepeat = 1:10
		disp(sprintf('%d/%d, %d/%d repeat\n', iN_hat, numel(N_hat_vals), irepeat, 10));

		n = sample_chain(rateFunc, serviceDistn, T, N_hat);
		y = sample_obs(alpha, n);
		disp(y)
		gamma = immigration_rate(rateFunc, serviceDistn, T, N_hat);
		delta = survival_prob(serviceDistn, T);

		start = tic;
		[a,b,c] = forward_messages(rateFunc, serviceDistn, T, y, N_hat, alpha, 'N_max', n_max);
		runtime_fa(irepeat, iN_hat) = toc(start);

		start = tic;
		[a,b,c,d,f,g] = gf_forward(y, gamma, alpha, delta);
		runtime_gffa(irepeat, iN_hat) = toc(start);
	end
end