function [output] = norm_proxg(norm, b,x, t, rho, gamma)
%DESCRIPTION: Prox operator for l1 and l2 norms 
%   INPUT:  norm             = problem (either l1 or l2)
%           b                = blurred image
%           t                = stepsize
%           rho              = relaxation parameter
%           x                = input vector
%           gamma            = de-noising parameter
%   OUTPUT: l1 or l2 depending on parameter
%
    if norm == l1
        %insert l1 prox g
    else
        %insert l2 prox g
    end
end