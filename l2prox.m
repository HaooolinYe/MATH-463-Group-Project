function [output] = l2prox(y1, b, t)
%DESCRIPTION: Prox operator of the iso norm 
%   INPUT:  
%           t                = stepsize
%           (y1,b)          = input vector
%   OUTPUT:  prox( l2(y1-b) )
%

x = 2*t*b+y1;
output = (x ./ (2*t + 1) );

end