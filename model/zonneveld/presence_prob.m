function p_t = presence_prob( arrivalDistn, serviceDistn, T )
% PRESENCE_PROB := Compute p_t, the probability an individual is alive at t for all t in T
% p_t = presence_prob( arrivalDistn, serviceDistn, T )
%
% INPUTS
% required:
%    arrivalDistn  = a distribution object for the birth process (typically normal)
%    serviceDistn  = a distribution object for the death process (typically exponential)
%    T             = vector [1 x K] of observation times (sample times)
%
% OUTPUTS
%    p_t           = vector [1 x K] of presence likelihood of an individual at all t in T
%                    note: this is _not_ at all a probability distribution

%for each t, compute p_t as the convolution over all possible birth times before t
%and the probability of living at least until t, computed from the cdf
p_t = arrayfun(@(t) ...
               quadgk(@(s) arrivalDistn.pdf(s) .* ...
               	           (1 - serviceDistn.cdf(t-s)), ...
                      -inf, t), ...
               T);

end

