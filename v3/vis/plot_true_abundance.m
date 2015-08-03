function plot_true_abundance( S, Z, varargin )

assert(size(S, 1) == size(Z, 1))

%defaults

parser = inputParser;
addOptional(  parser, 'T',         []);
addOptional(  parser, 'y',         []);

parser.parser(varargin{:})
T = parser.Results.T;
y = parser.Results.y;

%if S, Z are not cell arrays, wrap them in cell arrays
%using a 2D matrix would be claener, but would require every S,Z have the same # of indivs


end

