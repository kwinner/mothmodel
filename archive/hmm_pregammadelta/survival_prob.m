function delta = survival_prob( serviceDistn, T )
% SURVIVAL_PROB := evaluate the probability of an individual surviving from T(k) to T(k+1)
% delta = survival_prob( serviceDistn, T )
%
% INPUTS
% required:
%    serviceDistn = a distribution object for the death process (typically exponential)
%                   serviceDistn should be created with makedist(...)
%    T            = vector of observation times (reals)
%
% OUTPUTS
%    delta        = vector of survival probability of existing indivs

	%probability of an individual living from t1 to t2
	delta = arrayfun(@(t1,t2) (1-serviceDistn.cdf(t2 - t1)) ...
		             , T(1:end-1), T(2:end));

end

