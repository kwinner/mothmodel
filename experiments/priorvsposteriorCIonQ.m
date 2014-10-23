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

%do the prior CI
samples = zeros(MCMCsteps,T+1)
for isample = 1:MCMCsteps
	samples(isample,:) = binornd(params.N, Pb);
end

%compute the CI for the mean of a binomial with samples of each column
for iinterval = 1:T+1
end

%sample state from P
[y, state] = sampleState(theta, params);

%fit the prior distribution on B
%B(i) ~ Binomial(N, Pb(i));
[~, pci_prior] = binofit(B, params.N, alpha);

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
pci_post = B([1,end],:)';