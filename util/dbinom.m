%DBINOM Compute the binomial denisty (or log-density)
%
% y = dbinom(x, n, p, as_log, raw)
%
% Inputs:
%   x        vector of values at which to evalute the (log-) pmf
%   n        binomial parameter: number of trials
%   p        binomial parameter: success probability
%   as_log   if true, return log pmf                      (default: false)
%   raw      if true, don't force inputs to be integers   (default: false)
%
% Note that this is vectorized so that x may be a vector, but the 
% parameters n and p must be scalars,
%
% The implementation is based on C code from the lightspeed toolbox 
% by Tom Minka, and the Rmath library (with code by Catherine Loader).
% See comments in dbinom.c for full attribution and license.
% 