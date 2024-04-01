function [ output ] = proxl1(t, x)
%DESCRIPTION: Prox operator of the l1 norm 
%   INPUT:  
%           t                = stepsize
%           x                = input vector
%   OUTPUT: prox operator of l1 norm (soft thresholding) applied to x
    output = zeros(size(x));
    for k = 1:numel(x)
        if abs(x(k))>t
            output(k) = sign(x(k)) * (abs(x(k))-t);
        end
    end
end