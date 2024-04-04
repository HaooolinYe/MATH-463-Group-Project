function [ output ] = proxl1(t, b, y1)
%DESCRIPTION: Prox operator of the l1 norm 
%   INPUT:  
%           t                = stepsize
%           b                = blurred image
%           x                = input vector
%   OUTPUT: prox operator of l1 norm (soft thresholding) applied to x
    z = y1-b;
    logic = abs(z)>t;
    output = logic.*(sign(z).*(abs(z)-t))+b;
end