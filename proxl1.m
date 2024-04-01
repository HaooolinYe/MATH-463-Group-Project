function [ output ] = proxl1(t, x)
%DESCRIPTION: Prox operator of the l1 norm 
%   INPUT:  
%           t                = stepsize
%           x                = input vector
%   OUTPUT: prox operator of l1 norm (soft thresholding) applied to x
    absxk = abs(x);
    logic = absxk>t;
    output = logic.*(sign(x).*(absxk-t));
end