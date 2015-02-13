function NLL = observation_NLL( y, N, p_t, alpha )
%compute the log likelihood of the observations given the abundancy values

NLL = -sum(logpoisspdf(y, N .* p_t .* alpha));

end

