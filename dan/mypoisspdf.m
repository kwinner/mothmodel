% Simulate calculation of Poisson pdf using gf_forward
function p = mypoisspdf(lambda, y, k)

    % Poisson arrivals at first time step; no immigration
    gamma = [lambda zeros(1,k-1)];

    % 100% survival
    delta = ones(1,k-1);
    
    % 100% detection
    alpha = 1;

    % same observation of y at every time-step
    y = repmat(y, 1, k);

    p = gf_forward( y, gamma, alpha, delta );
end