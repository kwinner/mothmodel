N      = 100;
mu     = 5;
sigma  = 1;
lambda = 1;
alpha  = .5;
t      = 1:10;
T      = numel(t);
NQ     = (T+1)*(T+2)/2;

N_ITER = 1000;


ysample  = cell(1,N_ITER);
qsample  = cell(1,N_ITER);
nsample  = cell(1,N_ITER);
psample  = cell(1,N_ITER);
LLsample = cell(1,N_ITER);
qgibbs   = cell(1,N_ITER);
ngibbs   = cell(1,N_ITER);
LLgibbs  = cell(1,N_ITER);

for i = 1:N_ITER
    [ysample{i}, qsample{i}, nsample{i}, psample{i}] = generateMothData(N, mu, sigma, lambda, alpha, t);
    LLsample{i} = mothLikelihood(psample{i},qsample{i},nsample{i},ysample{i},alpha);
    
    DELTA = randsample(NQ,2,false);
    [a,b] = ind2sub_triu(T+1, DELTA(1));
    [c,d] = ind2sub_triu(T+1, DELTA(2));    
    [qgibbs{i}, ngibbs{i}] = gibbs(psample{i},qsample{i},nsample{i},ysample{i},alpha,a,b,c,d,'logspace',true);
    LLgibbs{i} = mothLikelihood(psample{i},qgibbs{i},ngibbs{i},ysample{i},alpha);
end