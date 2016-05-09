function [beta_gf, gamma_gf, psi_gf] = gf_backward( y, arrivRate, detectProb, deathProb )
% GF_BACKWARD := backward algorithm for count models implemented with generating functions
% beta = gf_backward( y, gamma, alpha, delta )
%
% INPUTS
% required:
%    y     = vector [1 x K] of observed counts
%    arrivRate  = vector [1 x K] of arrival rate
%    detectProb = detection probability
%    deathProb  = vector [1 x K] of death probability
%
% OUTPUTS
%    beta = generating functions for backward messages

K = length(arrivRate);

syms s ev;
beta_gf(K)   = symfun((1-s)^-1, s);
betaDeriv(K) = diff(beta_gf(K), y(K), s);
gamma_gf(K)  = symfun(1/(gamma(y(K)+1))*(s*detectProb)^(y(K)) * subs(betaDeriv(K), s, s*(1-detectProb)), s);

for k = K-1:-1:1
	arrivRate_k = arrivRate(k);
	deathProb_k = deathProb(k);
	psi_gf(k) = symfun(gamma_gf(k+1)*exp(arrivRate_k*(s^-1-1)),s);
	beta_gf(k) = symfun(1/(deathProb_k * s - s + 1) * subs(psi_gf(k), s, (deathProb_k*s/(deathProb_k*s-s+1))), s);

	betaDeriv(k) = diff(beta_gf(k), y(k), s);
	gamma_gf(k)  = symfun(1/(gamma(y(k)+1))*(s*detectProb)^(y(k)) * subs(betaDeriv(k), s, s*(1-detectProb)), s);
end

%evaluate the GF at 1
% likelihood = sum(f .* ((c + d) .^ (0:length(f)-1))) * exp(a + b);

end

