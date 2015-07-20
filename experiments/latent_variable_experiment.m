function latent_variable_experiment()

no_plot = false;

%parameters learned for population A1 in the Gross et al paper
mu     = 6.9;
sigma  = 3.4;
lambda = 1/0.645;
alpha  = .32;

if alpha == 1
    N = 70;
else
    N = 216;
end

%load data
data = sfs();
t = [data(~isnan([data.A1])).t];
y = [data(~isnan([data.A1])).A1];

%compute canonical state and initial P
Q_0 = canonical_state(y, N);
P = birthdeath_pmf(mu, sigma, lambda, N, t);

moves  = {'pair','shuffle','mergesplit'};
mcmc_nIter = 2000;
use_ARS = false;
num_runs = 2;

%whether to draw samples from prior or posterior
prior = false;

postfix = '_posterior';
if prior
    postfix = '_prior';
end

% if ~no_plot
%     figure; hold on;
% end
cum_vals_vec = zeros(num_runs,2*N);
times_vec = zeros(num_runs,2*N);
y_vec = zeros(num_runs,length(t));
S = zeros(N,num_runs);
Z = zeros(N,num_runs);
for n = 1:num_runs
    fprintf('%d\n', n);
    if prior
        Q = Q_0;
    else
        [~, Q_history] = mcmc(y, Q_0, P, N, alpha, t, 'nIterations',mcmc_nIter,'moves',moves,'useARS',logical(use_ARS));
        Q = Q_history{end};
    end
    
    draws = sample_from_Q(Q,mu,sigma,lambda,N,t);
    S(:,n) = draws(:,1);
    Z(:,n) = draws(:,2) - draws(:,1);

    %set up observed curve
    vals = [ones(N,1); -1*ones(N,1)];
    times = [draws(:,1); draws(:,2)];
    [times,inds] = sort(times);
    vals = vals(inds);
    cum_vals = cumsum(vals);
    
    cum_vals_vec(n,:) = cum_vals;
    times_vec(n,:) = times;
    y_vec(n,:) = y;
    
    % if ~no_plot
    %     stairs(times,cum_vals,'edgealpha',.2)
    % end
end

plot_true_abundance(S,Z,t,y, 'lineColor', [0, 0.447, 0.741, 1]);
plot(min(t):.1:max(t), N.*presence_probs(mu,sigma,lambda,min(t):.1:max(t)));
keyboard
return

%get pdf
ks = floor(min(times)):ceil(max(times));
ps = zeros(size(ks));
for k_ind = 1:length(ks)
    k = ks(k_ind);
    %     func = @(s) 1/(sigma*sqrt(2*pi))*exp(-(s-mu).^2/(2*sigma^2)).*(1 - exp(-1/lambda*(k-s)));
    func = @(s) normpdf(s,mu,sigma).*expcdf(k-s,lambda,'upper');
    ps(k_ind) = integral(func, -inf, k);
end
hold on;
if ~no_plot
    plot(ks,ps*N,'r','linewidth',4);
    title(sprintf('N = %d',N))
end

csvwrite([fileparts(which(mfilename)) '/latent_variable_cum_vals',postfix,'.csv'],cum_vals_vec)
csvwrite([fileparts(which(mfilename)) '/latent_variable_times',postfix,'.csv'],times_vec)
csvwrite([fileparts(which(mfilename)) '/latent_variable_ys',postfix,'.csv'],y_vec)
csvwrite([fileparts(which(mfilename)) '/latent_variable_mean_ps_ks',postfix,'.csv'],[ps;ks])

keyboard

end

function draws = sample_from_Q(Q,mu,sigma,lambda,N,t)

% how many of each q bin we've filled so far
counts = zeros(size(Q));

% list of birth and death times
draws = zeros(N,2); % [birth time, death time]
iDraw = 1;

%append -inf and inf to t
tWithBounds = [-inf, t, inf];

for iBirth = 1:size(Q,1)
    %create a truncated birth distribution which forces births to be in the right window
    smin = tWithBounds(iBirth);
    smax = tWithBounds(iBirth + 1);
    birthDist = truncate(makedist('Normal', 'mu', mu, 'sigma', sigma), smin, smax);

    for iDeath = iBirth:size(Q,2)
        %truncate the lifespan distribution pessimistically
        %min lifespan = mindeathtime - maxbirthtime
        zmin = tWithBounds(iDeath) - smax;
        %max lifespan = maxdeathtime - minbirthtime
        zmax = tWithBounds(iDeath + 1) - smin;

        lifespanDist = truncate(makedist('Exponential', 'mu', lambda), zmin, zmax);

        while counts(iBirth, iDeath) < Q(iBirth, iDeath)
            %draw a birth and death
            birth    = random(birthDist);
            lifespan = random(lifespanDist);
            death    = birth + lifespan;

            %check if the individual is needed in the counts
            iDeath = max(find(t < death));
            if counts(iBirth, iDeath) < Q(iBirth, iDeath)
                draws(iDraw, :) = [birth, death];
                counts(iBirth, iDeath) = counts(iBirth, iDeath) + 1;
                iDraw = iDraw + 1;
            end
        end
    end

    while any(counts(iBirth, :) < Q(iBirth, :))
        %draw a birth and death
        birth    = random(birthDist);
        lifespan = random(lifespanDist);
        death    = birth + lifespan;

        %check if the individual is needed in the counts
        iDeath = max(find(t < death));
        if counts(iBirth, iDeath) < Q(iBirth, iDeath)
            draws(iDraw, :) = [birth, death];
            counts(iBirth, iDeath) = counts(iBirth, iDeath) + 1;
            iDraw = iDraw + 1;
        end
    end
end

% while ~all(all(counts == Q)) % until we've filled up counts to match q
    
%     % sample birth and death time
%     birth = normrnd(mu,sigma);
%     death = birth + exprnd(lambda);
    
%     %determine birth and death intervals
%     birth_interval = min(find(t >= birth)) - 1;
%     death_interval = max(find(t <= death));
    
%     %increment count if its not full
%     if counts(birth_interval+1,death_interval+1) < Q(birth_interval+1,death_interval+1)
%         counts(birth_interval+1,death_interval+1) = counts(birth_interval+1,death_interval+1) + 1;
        
%         %add to draws
%         draws(draw_ind,:) = [birth,death];
%         draw_ind = draw_ind + 1;
%     end
    
end
