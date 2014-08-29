N      = 100;
mu     = 5;
sigma  = 1;
lambda = 1;
alpha  = .5;
t      = 1:10;
T      = numel(t);
NQ     = (T+1)*(T+2)/2;

N_ITER = 1000;

q = cell(1,N_ITER);
n = cell(1,N_ITER);
normalizedLikelihood1 = cell(1,N_ITER);
normalizedLikelihood2 = cell(1,N_ITER);
normalizedLikelihood3 = cell(1,N_ITER);

[y, q{1}, n{1}, p] = generateMothData(N, mu, sigma, lambda, alpha, t);

n{1} = y;
n{1}(1) = N-sum(y);
q{1} = zeros(T+1);
for i = 1:T
    q{1}(i,i+1) = n{1}(i);
end
    
i = 1;
LL = zeros(N_ITER,3);
LL(i,:) = mothLikelihood(p,q{1},n{1},y,alpha);

for i=2:N_ITER
    DELTA = randsample(NQ,2,false);
    [a,b] = ind2sub_triu(T+1, DELTA(1));
    [c,d] = ind2sub_triu(T+1, DELTA(2));
    
    [q{i}, n{i}, normalizedLikelihood1{i}] = gibbs(p,q{i-1},n{i-1},y,alpha,a,b,c,d);
    [q{i}, n{i}, normalizedLikelihood2{i}] = gibbs(p,q{i-1},n{i-1},y,alpha,a,b,c,d,'logspace',true);
    [q{i}, n{i}, normalizedLikelihood3{i}] = gibbs(p,q{i-1},n{i-1},y,alpha,a,b,c,d,'likelihoodMode','full','logspace',true);
    LL(i) = mothLikelihood(p,q{i},n{i},y,alpha);
end

plot(-LL)

sum(cellfun(@(x,y) sum(abs(x-y)), normalizedLikelihood1, normalizedLikelihood2))
sum(cellfun(@(x,y) sum(abs(x-y)), normalizedLikelihood2, normalizedLikelihood3))
sum(cellfun(@(x,y) sum(abs(x-y)), normalizedLikelihood1, normalizedLikelihood3))
max(cell2mat(cellfun(@(x,y) max(abs(x-y)), normalizedLikelihood1, normalizedLikelihood2, 'UniformOutput', false)))
max(cell2mat(cellfun(@(x,y) max(abs(x-y)), normalizedLikelihood2, normalizedLikelihood3, 'UniformOutput', false)))
max(cell2mat(cellfun(@(x,y) max(abs(x-y)), normalizedLikelihood1, normalizedLikelihood3, 'UniformOutput', false)))