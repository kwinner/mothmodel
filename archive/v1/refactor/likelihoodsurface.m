MU = 5:1:13;    mu_0 = 9.5;
SIGMA = [.1:.1:8];   sigma_0 = 2.8;
LAMBDA = [0.1:.25:15];  lambda_0 = 7;

alpha = 1;
N     = 50;

%realdata = sfs;
%t = [realdata(~isnan([realdata.E2])).t] + 1;
%y = [realdata(~isnan([realdata.E2])).E2];
t = 0:5:40;
y = sampleState(struct('mu', mu_0, 'sigma', sigma_0, 'lambda', lambda_0), struct('N', N, 'alpha', alpha, 't', t));

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
LLzonn_mulambda = zeros(numel(MU), numel(LAMBDA));
LLGP_mulambda = zeros(numel(MU), numel(LAMBDA));
cnt = 1;
tot = numel(LLzonn_mulambda);
for imu = 1:numel(MU)
	for ilambda = 1:numel(LAMBDA)
        tic
        mu     = MU(imu);
		lambda = LAMBDA(ilambda);
		N      = N_RANGE(iN);
		LLzonn_mulambda(imu, ilambda) = -zonn_loglikelihood(y, [mu, sigma_0], lambda, t, alpha, N);
		LLGP_mulambda(imu, ilambda)   = -GPLL(struct('mu',mu,'sigma',sigma_0,'lambda',lambda), y, struct('N',N,'alpha',alpha,'t',t));

        fprintf('Iter %d/%d done in %.2f seconds\n', cnt, tot, toc);
        cnt = cnt + 1;
    end
end

figure
hold on
surf(LAMBDA, N_RANGE, LLzonn_lambdaN');
plot3(lambda_0, N_0, -zonn_loglikelihood(y, [mu_0, sigma_0], lambda_0, t, alpha, N_0));
xlabel('Lambda');
ylabel('N');
title('zonn')

figure
hold on
surf(LAMBDA, N_RANGE, LLGP_lambdaN');
plot3(lambda_0, N_0, -GPLL(struct('mu',mu_0,'sigma',sigma_0,'lambda',lambda_0), y, struct('N',N_0,'alpha',alpha,'t',t)));
xlabel('Lambda');
ylabel('N');
title('GP')

hold off