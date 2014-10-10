function LL = zonneveldLL(theta, y, N, alpha, t)
mu = theta(1); sigma = theta(2); lambda = theta(3);

%compute the matrix p (a multinomial with compact parameters theta)
p = pdf_p(theta, t);

%for each observation, compute the probability of being alive at that time
%note, this is mathematically equivalent to finding abundance on p instead of q
T = size(p,1) - 1;

%from abundancy.mat, but real abundancy needs to cut out rounding error
%each element of n is a 2D sum of a submatrix of q
pt = arrayfun(@(k) sum(sum(p(1:k, k+1:T+1))), 1:T);

%convert pt from "prob alive at t_k" to "prob observed at t_k"
pt = alpha .* pt;

%compute the LL
LL =    logfactorial(N)   ...
      - logfactorial(y)   ...
      - logfactorial(N-y) ...
      + y .* log(pt)      ...
      + (N-y) .* log(1-pt);

%sum the LL of each observation
LL = sum(LL);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                           PDF for cell prob as a func of theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P] = pdf_p (theta, t)
%for integration, use -inf and inf as min and max vals on t
T = numel(t);
t = [-inf, t, inf];
P = zeros(T+1);
for i = 1:T+1
	for j = i:T+1
		P(i,j) = quadgk(@(s) pdf_p_integrand(i, j, s, theta, t), t(i), t(i+1));
	end
end
end
function [integrand] = pdf_p_integrand (i, j, s, theta, t)
mu = theta(1); sigma = theta(2); lambda = theta(3);
[zmin, zmax] = lifespan_domain(t(j), t(j+1), s);
integrand = normpdf(s, mu, sigma) .* (expcdf(zmax, lambda) - expcdf(zmin, lambda));
end
%% utility function
% there's a problem with the integration of F(zmax) - F(zmin)
% F(z) is 0 for z <= 0, even though z \in [0,inf)
% this isn't a problem for the pdf, but the gradients need to restrict z appropriately
% especially dp/dlambda, which includes multiplication by zmin/zmax
function [zmin, zmax] = lifespan_domain (death_min, death_max, birth)
zmin = death_min - birth;
zmin(death_min == birth) = 0;
zmin = max(zmin, 0); %clamp z to 0
zmax = death_max - birth;
zmax(death_max == birth) = inf;
end