function gamma = immigration_rate( rateFunc, serviceDistn, T, N_hat )
% IMMIGRATION_RATE := evaluate the rate of immigration between each observation time
% gamma = immigration_rate( rateFunc, serviceDistn, T )
%
% INPUTS
% required:
%    rateFunc     = the arrival rate function for the time-varying Poisson arrival process (typically normal)
%    serviceDistn = a distribution object for the death process (typically exponential)
%                   serviceDistn should be created with makedist(...)
%    T            = vector of observation times (reals)
%    N_hat        = superpopulation size
%
% OUTPUTS
%    gamma        = vector of rates of new "successful" immigrants

%rate of arrivals who survive to t_0
gamma = N_hat .* arrayfun(@(t1,t2) ...
	                               quadgk(@(s) ...
	                                      rateFunc(s) ...
	                                      .* (1-serviceDistn.cdf(t2 - s)) ...
	                                      , t1, t2) ...
	                      , [-inf, T(1:end-1)], T);

end

