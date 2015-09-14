function latent_variable_experiment()

no_plot = true;

%parameters learned for population A1 in the Gross et al paper
mu     = 6.9;
sigma  = 3.4;
lambda = 1/0.645;
alpha  = 1;%.32;
N      = 70;
t      = 1:20;
params = struct('N', N, 'alpha', alpha, 't', t);
theta = struct('mu', mu, 'sigma', sigma, 'lambda', lambda);

moves  = {'pairwise','shuffle','mergesplit'};
mcmc_nIter = 1000;
use_ARS = false;
num_runs = 3;

%whether to draw samples from prior or posterior
prior = false;

postfix = '_posterior';
if prior
    postfix = '_prior';
end

if ~no_plot
    figure; hold
end
cum_vals_vec = zeros(num_runs,2*N);
times_vec = zeros(num_runs,2*N);
y_vec = zeros(num_runs,length(t));
for n = 1:num_runs
    
    % get q
    [y, state_0] = sampleState(theta,params);
    
    if prior == true
        q = state_0.q;
    else
        state = mcmc(y,state_0,params,'iterations',mcmc_nIter,'moves',moves,'use_ARS',logical(use_ARS));
        q = state(end).q;
    end
    
    draws = sample_from_q(q,params,theta);
    
    %set up observed curve
    vals = [ones(N,1); -1*ones(N,1)];
    times = [draws(:,1); draws(:,2)];
    [times,inds] = sort(times);
    vals = vals(inds);
    cum_vals = cumsum(vals);
    
    cum_vals_vec(n,:) = cum_vals;
    times_vec(n,:) = times;
    y_vec(n,:) = y;
    
    if ~no_plot
        stairs(times,cum_vals,'edgealpha',.2)
    end
end

%get pdf
ks = floor(min(times)):ceil(max(times));
ps = zeros(size(ks));
for k_ind = 1:length(ks)
    k = ks(k_ind);
    %     func = @(s) 1/(sigma*sqrt(2*pi))*exp(-(s-mu).^2/(2*sigma^2)).*(1 - exp(-1/lambda*(k-s)));
    func = @(s) normpdf(s,mu,sigma).*expcdf(k-s,lambda,'upper');
    ps(k_ind) = integral(func, -inf, k);
end

if ~no_plot
    plot(ks,ps*N,'linewidth',4);
    title(sprintf('N = %d',N))
end

csvwrite([fileparts(which(mfilename)) '/latent_variable_cum_vals',postfix,'.csv'],cum_vals_vec)
csvwrite([fileparts(which(mfilename)) '/latent_variable_times',postfix,'.csv'],times_vec)
csvwrite([fileparts(which(mfilename)) '/latent_variable_ys',postfix,'.csv'],y_vec)
csvwrite([fileparts(which(mfilename)) '/latent_variable_mean_ps_ks',postfix,'.csv'],[ps;ks])



function draws = sample_from_q(q,params,theta)

t = params.t;
N = params.N;
mu = theta.mu;
sigma = theta.sigma;
lambda = theta.lambda;

% how many of each q bin we've filled so far
counts = zeros(size(q));

% list of birth and death times
draws = zeros(N,2); % [birth time, death time]
draw_ind = 1;

while ~all(all(counts == q)) % until we've filled up counts to match q
    
    % sample birth and death time
    birth = normrnd(mu,sigma);
    death = birth + exprnd(lambda);
    
    %determine birth and death intervals
    birth_interval = min(max(floor(birth), 0), max(t));
    death_interval = min(max(floor(death), 0), max(t));
    
    %increment count if its not full
    if counts(birth_interval+1,death_interval+1) < q(birth_interval+1,death_interval+1)
        counts(birth_interval+1,death_interval+1) = counts(birth_interval+1,death_interval+1) + 1;
        
        %add to draws
        draws(draw_ind,:) = [birth,death];
        draw_ind = draw_ind + 1;
    end
    
end
