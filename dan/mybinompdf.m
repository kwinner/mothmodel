% Simulate calculation of Poisson pdf using gf_forward
function prob = mybinompdf(k, n, p)

    % Poisson arrivals at first time step
    gamma = [n 0];

    % survival p
    delta = p;
    
    % 100% detection
    alpha = 1;

    % Observations n, k
    y = [n k];

    prob = gf_forward( y, gamma, alpha, delta ) / poisspdf(n, n);
end