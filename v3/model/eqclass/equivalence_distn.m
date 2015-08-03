function P = equivalence_distn( arrivalDistn, serviceDistn, T )
% EQUIVALENCE_DISTN := Compute P, the multinomial distribution for Q, where P(i,j) = p(born in interval i, die in interval j)
% P = equivalence_distn( arrivalDistn, serviceDistn, T )
%    arrivalDistn  = a distribution object for the birth process (typically normal)
%    serviceDistn  = a distribution object for the death process (typically exponential)
%    T             = vector [1 x K] of observation times (sample times)
%
%    P             = probability matrix [K+1 x K+1] for likelihood of each equiv class
%                    organized as: rows = birth interval, cols = death interval

%readability definitions
K = size(T, 2);

%initialize outputs
P = zeros(K+1);

%convert sampling times, T, to interval bounds, which cover (-inf, inf)
intBounds = [-inf, T, inf];

%loop over the upper triangular part of P
stabilityWarning = false; %do we need to warn the user of a stability correction?
for i = 1:K+1
	for j = i:K+1
		%compute the convolution of arrival and lifespan processes for this cell
		%see section 3 of our ICML paper (Winner et al, 2015) for more details
		P(i,j) = quadgk(@(s) arrivalDistn.pdf(s) .* ...
			                 (serviceDistn.cdf(intBounds(j+1) - s) - ...
			                  serviceDistn.cdf(intBounds(j) - s)), ...
			            intBounds(i), intBounds(i+1));

		%none of the upper triangular outcomes should actually have probability 0
		%any zeros are the result of numerical inaccuracy and will lead to instability,
		%so we give them the minimum probability possible
		if P(i,j) == 0
			P(i,j) = realmin;
			stabilityWarning = true;
		end
	end
end

if stabilityWarning
	warning('Some elements of P set to realmin, check correctness of P before using...')
end

end

