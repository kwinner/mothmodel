function [p, p_n_y ] = dail_madsen_trunc(lambda, alpha, y, N)

    n = (0:N)';      % column vector
    y = y(:)';       % row vector
    k = length(y);
    
    p_n = poisspdf(n, lambda); % row vector
    p_y_given_n = binopdf( repmat(y, N+1, 1), repmat(n, 1, k), alpha); % N x k matrix

    p_n_y = p_n .* prod(p_y_given_n, 2);
    
    p = sum( p_n_y, 1);

end
