function [ error ] = relativeerror( theta, theta_hat )
%RMSE Summary of this function goes here
%   Detailed explanation goes here

if isstruct(theta)
	theta = [theta.mu, theta.sigma, theta.lambda];
end
if isstruct(theta_hat)
	theta_hat = [theta_hat.mu, theta_hat.sigma, theta_hat.lambda];
end

error = sqrt(mean((theta_hat - theta).^2./theta));

end

