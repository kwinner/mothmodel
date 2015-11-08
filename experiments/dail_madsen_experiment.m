function dail_madsen_experiment()

% runtime_ll_vs_nmax_experiment();
% runtime_vs_alpha_experiment();
% runtime_vs_y_experiment();
% runtime_vs_K_experiment();
runtime_vs_gamma_experiment();

end

function runtime_vs_gamma_experiment()
gamma = [100:100:2000];
K = 10;
alpha = .4;

nrepeats = 20;

n = poissrnd(1000);
y = binornd(repmat(n,1,K), alpha);

gffa_runtime = zeros(1, length(gamma));
dm_runtime = zeros(1, length(gamma));
for igamma = 1:length(gamma)
	disp(gamma(igamma))
	N_max = 1.25 * gamma(igamma);
	for irepeat = 1:nrepeats
		tic
		gffa_likelihood = gf_forward(y, [gamma(igamma),zeros(1,K-1)], alpha, 1);
		gffa_runtime(igamma) = gffa_runtime(igamma) + toc;
		tic
		dm_likelihood = dail_madsen(y, gamma(igamma), alpha, N_max);
		dm_runtime(igamma) = dm_runtime(igamma) + toc;
	end
end

gffa_runtime = gffa_runtime ./ nrepeats;
dm_runtime   = dm_runtime   ./ nrepeats;

figure
plot(gamma, dm_runtime, gamma, gffa_runtime, 'LineWidth', 4)
title('Runtime vs gamma','Fontsize',20)
xlabel('gamma','FontSize',16)
ylabel('Runtime (s)', 'FontSize',16)
xlim([min(gamma), max(gamma)])
legend({'Dail-madsen', 'GF-FA'},'FontSize',14,'Location','northeast')

end

function runtime_vs_K_experiment()
gamma = 100;
K = 1:25;
alpha = .6;

nrepeats = 20;

N_max = 1.25 * gamma;

n = poissrnd(gamma);

gffa_runtime = zeros(1, length(K));
dm_runtime = zeros(1, length(K));
for iK = 1:length(K)
	disp(K(iK))
	y = binornd(repmat(n,1,K(iK)), alpha);
	for irepeat = 1:nrepeats
		tic
		gffa_likelihood = gf_forward(y, [gamma,zeros(1,K(iK)-1)], alpha, 1);
		gffa_runtime(iK) = gffa_runtime(iK) + toc;
		tic
		dm_likelihood = dail_madsen(y, gamma, alpha, N_max);
		dm_runtime(iK) = dm_runtime(iK) + toc;
	end
end

gffa_runtime = gffa_runtime ./ nrepeats;
dm_runtime   = dm_runtime   ./ nrepeats;

figure
plot(K, dm_runtime, K, gffa_runtime, 'LineWidth', 4)
title('Runtime vs K','Fontsize',20)
xlabel('K','FontSize',16)
ylabel('Runtime (s)', 'FontSize',16)
xlim([min(K), max(K)])
legend({'Dail-madsen', 'GF-FA'},'FontSize',14,'Location','northeast')

end

function runtime_vs_y_experiment()
gamma = 100;
K = 10;
alpha = .6;

nrepeats = 20;

N_max = 1.25 * gamma;
y = repmat(1:gamma,K,1)';

gffa_runtime = zeros(1, size(y,1));
dm_runtime = zeros(1, size(y,1));
for iy = 1:size(y,1)
	disp(y(iy,1))
	for irepeat = 1:nrepeats
		tic
		gffa_likelihood = gf_forward(y(iy,:), [gamma,zeros(1,K-1)], alpha, 1);
		gffa_runtime(iy) = gffa_runtime(iy) + toc;
		tic
		dm_likelihood = dail_madsen(y(iy,:), gamma, alpha, N_max);
		dm_runtime(iy) = dm_runtime(iy) + toc;
	end
end

gffa_runtime = gffa_runtime ./ nrepeats;
dm_runtime   = dm_runtime   ./ nrepeats;

figure
plot(y(:,1), dm_runtime, y(:,1), gffa_runtime, 'LineWidth', 4)
title('Runtime vs y','Fontsize',20)
xlabel('y','FontSize',16)
ylabel('Runtime (s)', 'FontSize',16)
xlim([min(y(:,1)), max(y(:,1))])
legend({'Dail-madsen', 'GF-FA'},'FontSize',14,'Location','northeast')

end

function runtime_vs_alpha_experiment()
gamma = 100;
K = 10;
alpha = 0.01:.01:1;

nrepeats = 10;

n = poissrnd(gamma);

N_max = 1.25 * gamma;

gffa_runtime = zeros(1, length(alpha));
dm_runtime = zeros(1, length(alpha));
for ialpha = 1:length(alpha)
	disp(alpha(ialpha))
	for irepeat = 1:nrepeats
		y = binornd(repmat(n,1,K), alpha(ialpha));
		tic
		gffa_likelihood = gf_forward(y, [gamma,zeros(1,K-1)], alpha(ialpha), 1);
		gffa_runtime(ialpha) = gffa_runtime(ialpha) + toc;
		tic
		dm_likelihood = dail_madsen(y, gamma, alpha(ialpha), N_max);
		dm_runtime(ialpha) = dm_runtime(ialpha) + toc;
		y = binornd(repmat(n,1,K), alpha(ialpha));
		tic
		dm_likelihood = dail_madsen(y, gamma, alpha(ialpha), N_max);
		dm_runtime(ialpha) = dm_runtime(ialpha) + toc;
		tic
		gffa_likelihood = gf_forward(y, [gamma,zeros(1,K-1)], alpha(ialpha), 1);
		gffa_runtime(ialpha) = gffa_runtime(ialpha) + toc;
	end
end

gffa_runtime = gffa_runtime ./ nrepeats;
dm_runtime   = dm_runtime   ./ nrepeats;

figure
plot(alpha, dm_runtime, alpha, gffa_runtime, 'LineWidth', 4)
title('Runtime vs \alpha','Fontsize',20)
xlabel('\alpha','FontSize',16)
ylabel('Runtime (s)', 'FontSize',16)
xlim([min(alpha), max(alpha)])
legend({'Dail-madsen', 'GF-FA'},'FontSize',14,'Location','northeast')

end

function runtime_ll_vs_nmax_experiment()
%standard
gamma = 100;
K = 10;
alpha = 0.35;
nrepeats = 20;

%fast
% gamma = 10;
% K = 3;
% alpha = 0.4;
% nrepeats = 1;

n = poissrnd(gamma);
y = binornd(repmat(n,1,K), alpha);

N_max = max(y):1:1.5*gamma;
% N_max = [1.25*gamma,1.5*gamma];

gffa_runtime = 0;
for irepeat = 1:nrepeats
	tic
	gffa_likelihood = gf_forward(y, [gamma,zeros(1,K-1)], alpha, 1);
	gffa_runtime = gffa_runtime + toc;
end

dm_runtime = zeros(1,length(N_max));
for i_nmax = 1:length(N_max)
	disp(N_max(i_nmax))
	for irepeat = 1:nrepeats
		tic
		dm_likelihood(i_nmax) = dail_madsen(y, gamma, alpha, N_max(i_nmax));
		dm_runtime(i_nmax) = dm_runtime(i_nmax) + toc;
	end
end

gffa_runtime = gffa_runtime ./ nrepeats;
dm_runtime   = dm_runtime   ./ nrepeats;

figure
plot(N_max, dm_runtime, [min(N_max), max(N_max)], [gffa_runtime, gffa_runtime], 'LineWidth', 4)
title('Runtime vs Problem Size','Fontsize',20)
xlabel('N_{max}','FontSize',16)
ylabel('Runtime (s)', 'FontSize',16)
xlim([min(N_max), max(N_max)])
legend({'Dail-madsen', 'GF-FA'},'FontSize',14,'Location','northwest')

figure
plot(N_max, dm_likelihood, [min(N_max), max(N_max)], [gffa_likelihood, gffa_likelihood], 'LineWidth', 4)
title('Likelihood Approximation','Fontsize',20)
xlabel('N_{max}','FontSize',16)
ylabel('Likelihood', 'FontSize',16)
xlim([min(N_max), max(N_max)])
legend({'Dail-madsen', 'GF-FA'},'FontSize',14,'Location','southeast')

end