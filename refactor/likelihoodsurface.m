MU = [-200:1:50];  mu_0 = 6.9;
SIGMA = [.1:.1:8];   sigma_0 = 3.4;
LAMBDA = [1:1:200]; lambda_0 = 1/.645;

alpha = .32;
N     = 216;

realdata = sfs;
t = [realdata(~isnan([realdata.A1])).t] + 1;
y = [realdata(~isnan([realdata.A1])).A1];

%hold mu fixed, vary sigma, lambda
% LL_sigmalambda = zeros(numel(SIGMA), numel(LAMBDA));
% for isigma = 1:numel(SIGMA)
% 	for ilambda = 1:numel(LAMBDA)
% 		sigma  = SIGMA(isigma);
% 		lambda = LAMBDA(ilambda);

% 		LL_sigmalambda(isigma, ilambda) = -zonn_loglikelihood(y, [mu_0, sigma], lambda, t, alpha, N);
% 	end
% end

% surf(LL_sigmalambda);
% xlim([min(SIGMA) max(SIGMA)]);
% xlabel('Sigma');
% ylim([min(LAMBDA) max(LAMBDA)]);
% ylabel('Lambda');

%hold sigma fixed, vary mu, lambda
LL_mulambda = zeros(numel(MU), numel(LAMBDA));
for imu = 1:numel(MU)
	for ilambda = 1:numel(LAMBDA)
		mu     = MU(imu);
		lambda = LAMBDA(ilambda);

		LL_mulambda(imu, ilambda) = -zonn_loglikelihood(y, [mu, sigma_0], lambda, t, alpha, N);
	end
end

surf(MU, LAMBDA, LL_mulambda');
xlabel('Mu');
ylabel('Lambda');