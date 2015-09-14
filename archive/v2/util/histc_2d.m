function [bincounts, ind_x, ind_y] = histc_2d( X, Y, x_edges, y_edges )

K = numel(x_edges);
N = numel(X);

[x_mesh, y_mesh] = meshgrid(1:K-1, 1:K-1);

%put all the instances into their bins
bins = arrayfun(@(i,j) find(X >= x_edges(i) & X < x_edges(i+1) & ...
                            Y >= y_edges(j) & Y < y_edges(j+1)), ...
                x_mesh', y_mesh', 'UniformOutput', false);

%count the number of elements in each bin
bincounts = cellfun(@numel, bins);

%optionally, report which bin each individual fell into
if nargout > 1
	ind_x = arrayfun(@(x) find(x >= x_edges & [x < x_edges(2:end),0]), X);
	ind_y = arrayfun(@(y) find(y >= y_edges & [y < y_edges(2:end),0]), Y);
end

end

