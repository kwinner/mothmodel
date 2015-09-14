mu = 9.5;
sigma = 2.8;
lambda = 7.0;
T_min = 0;
T_max = 50;


Ns = [25, 50, 100, 200];
T_freqs = [2, 5, 10];
alphas = [0.5, 1.0];
% Ns = [25, 50];
% T_freq = 10;
% alphas = 1.0;

for N = Ns
	for T_freq = T_freqs
		for alpha = alphas
			fprintf('N%d_u%d_a%.2f\n', N, T_freq, alpha)
			coverage_experiment(T_min, T_max, T_freq, mu, sigma, lambda, N, alpha, 100, sprintf('N%d_u%d_a%.2f', N, T_freq, alpha));
		end
	end
end