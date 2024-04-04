function [ output ] = proxl1(t, b, y1)
%DESCRIPTION: Prox operator of the l1 norm 
%   INPUT:  
%           t                = stepsize
%           b                = blurred image
%           x                = input vector
%   OUTPUT: prox operator of l1 norm (soft thresholding) applied to x
    z = b - y1;
    absxk = abs(z);
    logic = absxk>t;
    output = logic.*(sign(z)*t+y1);
end