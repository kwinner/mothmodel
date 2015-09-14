function [binCounts, indX, indY] = histc_2d( X, Y, xEdges, yEdges )
% HISTC_2D := Compute bin counts in two dimensions
% binCounts               = histc_2d( X, Y, xEdges, yEdges )
%                           (compute just the binCounts)
% [binCounts, indX, indY] = histc_2d( X, Y, xEdges, yEdges )
%                           (also compute the bin membership of each individual)
%                           note: at present, this is nontrivial
% INPUTS
% required:
%    X         = vector [N x 1] of data
%    Y         = vector [N x 1] of data
%    xEdges    = vector [1 x A] of bin bounds
%    yEdges    = vector [1 x B] of bin bounds
%          note: bin bounds should provide lower and upper bounds,
%                so append [-inf, *Edges, inf] to cover the full range
%          note: x bounds of bin {i,j} are [xEdges(i), xEdges(i+1))
%                and equivalently for y bounds
%                i.e. lower bounds are inclusive, upper bounds are not
%                (this pattern can be inverted by negating all input and transposing binCounts)
%
% OUTPUTS
%    binCounts = matrix [A-1 x B-1] of the number of instances which fell in each 2d bin
%    indX      = vector [N x 1] for each instance, the X index of its 2d bin
%    indY      = vector [N x 1] for each instance, the Y index of its 2d bin
%
% Original author unknown, adaptation by Kevin Winner
% kwinner@cs.umass.edu

assert(size(X, 1) == size(Y, 1))

%compute the lower (*LB) and upper (*UB) bounds of each bin, on a grid
[xLB, yLB] = meshgrid(xEdges(1:end-1), yEdges(1:end-1));
[xUB, yUB] = meshgrid(xEdges(2:end),   yEdges(2:end));

%find the instances belonging to each bin
membership = arrayfun(@(xlb,xub,ylb,yub) ...
	                  find(X >= xlb & X < xub & ...
                           Y >= ylb & Y < yub), ...
                      xLB, xUB, yLB, yUB, 'UniformOutput', false);

%count the number of instances in each bin
binCounts = cellfun(@numel, membership);

%optionally, report which bin each individual fell into
if nargout > 1
	ind_x = arrayfun(@(x) find(x >= xLB & x < xUB), X);
	ind_y = arrayfun(@(y) find(y >= yLB & y < yUB, Y);
end

end

