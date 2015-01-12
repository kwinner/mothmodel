function [ data ] = sfs ( filename )
%SFS Summary of this function goes here
%   Detailed explanation goes here

if ~exist('filename', 'var')
	filename = '/Users/kwinn/Work/Data/moth/sfs2004.csv';
end

data = csv2struct(filename, '%s%f%f%f%f');

%standardize the date
N = numel(data);
for i = 1:N
	%damnit this really should be doable with an arrayfun/structfun, come on matlab
	data(i).t = datenum(data(i).Date);
end

%diff by t_0
t_0 = min([data.t]);
for i = 1:N
	data(i).t = data(i).t - t_0;
end

data = rmfield(data, 'Date');

end

