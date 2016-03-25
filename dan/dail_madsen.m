function [p, a, b, c, d, f] = dail_madsen(lambda, alpha, y)

    k = length(y);
    
    % Poisson arrivals at first time step; no immigration
    gamma = [lambda zeros(1,k-1)];

    % 100% survival
    delta = ones(1,k-1);

    [p, a, b, c, d, f] = gf_forward( y, gamma, alpha, delta );

end
