function LL = latentvar_obs_LL( y, n, alpha, obs_dist )
%compute the log likelihood of the latent variables

if ~exist('obs_dist', 'var')
	obs_dist = 'bino';
end

switch obs_dist
case 'bino'
	LL = sum(dbinom(y, n, alpha, true));
case 'poiss'
	LL = sum(logpoisspdf(y, n .* alpha));
end


end

