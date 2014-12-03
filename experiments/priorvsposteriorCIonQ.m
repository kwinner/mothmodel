theta = struct( 'mu',     8,   ...
	            'sigma',  4,   ...
	            'lambda', 3    );
params = struct('N',      100, ...
	            'alpha',  1,   ...
	            't',      1:20 );
MCMCsteps    = 10000;
MCMCburnin   = 1000;
alpha        = .05;

T = numel(params.t);

%let B(i) be the total # of individuals born in interval i
%we can compute the CI on B(i) directly from the matrix P
%as follows:

P  = ppdf(theta, params); %matrix of all birth/death outcome probs
Pb = sum(P,1);            %vector of all birth probs

%sample state from P
[y, state] = sampleState(theta, params);

%fit the prior distribution on B
%B(i) ~ Binomial(N, Pb(i));
B_0 = sum(state.q,1);
[p_hat_prior, p_ci_prior] = binofit(B_0, params.N, alpha);
%convert from a CI on p to a CI on B
B_ci_prior = params.N .* p_ci_prior;

%compute a "CI" from the inverse CDF
B_ci_icdf_LB = binoinv(alpha/2, params.N, Pb);
B_ci_icdf_UB = binoinv(1-(alpha/2), params.N, Pb);
B_ci_icdf = [B_ci_icdf_LB ; B_ci_icdf_UB]';

%run mcmc
naive_state = naiveStateExplanation(y, params.N);
naive_state.p = P;
samples = mcmc(y, naive_state, params, 'iterations', MCMCsteps + MCMCburnin);

%remove burnin
samples = samples(MCMCburnin+1+1:end);

%compute B for each sample
B = zeros(MCMCsteps, T+1);
for isample = 1:MCMCsteps
	B(isample,:) = sum(samples(isample).q, 1);
end

%sort all the samples of B
B = sort(B,1);

%keep a CI worth of the B samples
CIsize = round(MCMCsteps * (1-alpha));
B = B(round((MCMCsteps - CIsize)/2):end-round((MCMCsteps - CIsize)/2),:);
B_ci_post = B([1,end],:)';

%okay now we have 3 CIs, let's compare their width
B_ci_prior_width = B_ci_prior(:,2) - B_ci_prior(:,1);
B_ci_icdf_width  = B_ci_icdf(:,2)  - B_ci_icdf(:,1);
B_ci_post_width  = B_ci_post(:,2)  - B_ci_post(:,1);

figure
plot(1:T+1, B_ci_icdf_width,  '-', ...
	 1:T+1, B_ci_post_width,  '-')
legend({'icdf','post'})
title('Width of confidence interval')
ylabel('CI width')
xlabel('Interval')

figure
hold on
h2 = plot(1:T+1, B_ci_icdf,  'color', 'green'); h2 = h2(1);
h3 = plot(1:T+1, B_ci_post,  'color', 'red');   h3 = h3(1);
h4 = plot(1:T+1, B_0, 'o',   'color', 'black')
legend([h2 h3],{'icdf','post'})
title('Births per interval')
ylabel('B')
xlabel('Interval')