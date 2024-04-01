function [output] = proxf(t, x)
%DESCRIPTION: Prox operator of f (indicator of S = {0<=x<=1})
%   INPUT:  t                = stepsize
%                   unused, kept for compatibility with other functions
%           x                = input vector
%   OUTPUT: prox operator of tf (projection onto S) applied to x
%
    output = max(min(1, x), 0); %projects each element in matrix to the set {0<=x<=1}
end


