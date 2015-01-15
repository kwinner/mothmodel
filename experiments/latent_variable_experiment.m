

%parameters learned for population A1 in the Gross et al paper
mu     = 6.9;
sigma  = 3.4;
lambda = 1/0.645;
alpha  = .32;
N      = 1000;
t      = 1:20;
params = struct('N', N, 'alpha', alpha, 't', t);
theta = struct('mu', mu, 'sigma', sigma, 'lambda', lambda);

% get q
[y, state] = sampleState(theta,params);
q = state.q;

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

%set up observed curve
vals = [ones(N,1); -1*ones(N,1)];
times = [draws(:,1); draws(:,2)];
[times,inds] = sort(times);
vals = vals(inds);
cum_vals = cumsum(vals);

%get pdf
ks = floor(min(times)):ceil(max(times));
ps = zeros(size(ks));
for k_ind = 1:length(ks)
    k = ks(k_ind);
%     func = @(s) 1/(sigma*sqrt(2*pi))*exp(-(s-mu).^2/(2*sigma^2)).*(1 - exp(-1/lambda*(k-s)));
    func = @(s) normpdf(s,mu,sigma).*expcdf(k-s,lambda,'upper');
    ps(k_ind) = integral(func, -inf, k);
end

figure; hold
stairs(times,cum_vals)
plot(ks,ps*N);
title(sprintf('N = %d',N))


