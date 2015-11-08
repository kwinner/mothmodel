function dan_test()

% Test 1: simulate Poisson pmf
lambda = 10;  % Poisson parameter
k = 10;       % 

y = (0:30)';
p1 = nan(size(y));
for i = 1:numel(y)
    p1(i) = mypoisspdf(lambda, y(i), k);
end
p2 = poisspdf(y, lambda);

plot(y, [p1 p2]);
max(abs(p1 - p2))


end

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