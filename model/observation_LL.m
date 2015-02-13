function LL = observation_LL( y, N, p_t, alpha )
%compute the log likelihood of the observations given the abundancy values

LL = sum(logpoisspdf(y, N .* p_t .* alpha));

end

