
%% Binomial
% n = 10;
% binomProb = 0.3;

% y = (0:10)';
% bn_p1 = nan(size(y));
% for i = 1:numel(y)
%     bn_p1(i) = mybinompdf(y(i), n, binomProb);
% end
% bn_p2 = binopdf(y, n, binomProb);

% plot(y, [bn_p1 bn_p2]);
% max(abs(bn_p1 - bn_p2))

% return

%% Dail Madsen
lambda = 1000;
alpha = 0.1;
y = [99 100 80];

tic
dm_p1 = dail_madsen(lambda, alpha, y)
toc

tic
dm_p2 = dail_madsen_trunc(lambda, alpha, y, 10000)
toc

return

%% Poisson pdf
% Test 1: simulate Poisson pmf
lambda = 10;  % Poisson parameter
k = 10;       % 

y = (0:30)';
gf_p1 = nan(size(y));
for i = 1:numel(y)
    gf_p1(i) = mypoisspdf(lambda, y(i), k);
end
gf_p2 = poisspdf(y, lambda);

plot(y, [gf_p1 gf_p2]);
max(abs(gf_p1 - gf_p2))

