function NLL = latentvar_obs_NLL( y, T, Q, alpha, obs_dist )
%compute the NLL of the latent variables

if ~exist('obs_dist', 'var')
	obs_dist = 'bino';
end

n = abundance(Q);

switch obs_dist
case 'bino'
	NLL = -sum(dbinom(y, n, alpha, true));
case 'poiss'
	NLL = -sum(logpoisspdf(y, n .* alpha));
end


end

