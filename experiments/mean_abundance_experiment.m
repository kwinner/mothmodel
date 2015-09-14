%initialize result storage
n_hmm = zeros(nSamples, numel(T));

%generate many samples from HMM
tStart = tic;
for iSample = 1:nSamples
	n_hmm(iSample, :) = sample_chain(rateFunc, serviceDistn, T, N_hat);
end
tEnd = toc(tStart);
runtimeHMM = tEnd;

%average the abundance at all observed times
n_mean_hmm = mean(n_hmm, 1);
n_std_hmm  = std(n_hmm, 1);

%compute mean abundance w/ Dan's code
tStart = tic;
n_mean_dan = N_hat .* dan_abundance(T, makedist('Normal','mu',mu,'sigma',sigma), {mu, sigma}, makedist('Exp','mu',lambda), {lambda});
tEnd = toc(tStart);
runtimeDan = tEnd;

figure; hold on
plot(T, n_mean_dan, '-', 'LineWidth', 3)
errorbar(T, n_mean_hmm, 2.*n_std_hmm, '--', 'LineWidth', 3)

%label it
title('Empirical mean abundance vs expected', 'FontSize', 20)
xlabel('time', 'FontSize', 16)
ylabel('abundance', 'FontSize', 16)
legend({'Expected','Empirical'})

% ylim([0, 1.05.*max(max(n_mean_hmm+2.*n_std_hmm),max(n_mean_dan))])