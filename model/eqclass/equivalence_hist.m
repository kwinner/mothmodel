function [Q, eqClass] = equivalence_hist( S, Z, T )
% EQUIVALENCE_HIST := Group individuals by equivalence class, then record counts of each class
% Q            = equivalence_hist( S, Z, T )
%                (compute the histogram only)
% [Q, eqClass] = equivalence_hist( S, Z, T )
%                (also compute the class membership of each individual)
%                note: this is non trivial with histc_2d and should be avoided
%
% INPUTS
% required:
%    S       = vector [N x 1] of individual birth times
%    Z       = vector [N x 1] of individual death times
%    T       = vector [1 x K] of observation times (sample times)
%
% OUTPUTS
%    Q       = matrix [K+1 x K+1] of counts for each equivalence class
%              organized as: rows = birth interval, cols = death interval
%    eqClass = vector [N x 1] of individual equivalence class memberships

assert(size(S, 1) == size(Z, 1));

%readability definitions
N = size(S, 1);
K = size(T, 2);
D = S + Z; %D = individual death times

if nargout == 1
	%avoid expensive index computation if eqClass unused
	Q = histc_2d(S, D, [-inf, T, inf], [-inf, T, inf]);
else
	[Q, ind_x, ind_y] = histc_2d(S, D, [-inf, T, inf], [-inf, T, inf]);

	%convert the x,y indices into equivalence classes for each individual
	[~,~,eqClass] = unique(sub2ind([K+1, K+1], ind_x, ind_y));
end

end

