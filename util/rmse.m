function [ error ] = rmse( y_hat, y )
%RMSE Summary of this function goes here
%   Detailed explanation goes here

if isstruct(y_hat)
    y_hat = struct2cell(y_hat);
    y_hat = cell2mat(y_hat);
    y_hat = mean(y_hat,3);
end

if isstruct(y)
    y = struct2cell(y);
    y = cell2mat(y);
    y = mean(y,3);
end

error = sqrt(mean((y - y_hat).^2));

end

