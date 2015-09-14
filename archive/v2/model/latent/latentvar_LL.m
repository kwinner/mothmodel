function LL = latentvar_LL( Q, P )
%compute the log likelihood of the latent variables

LL = logmnpdf(Q(:), P(:));

end

