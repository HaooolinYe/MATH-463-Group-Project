function [output] = iso_proxg(y2,y3, t, gamma)
%DESCRIPTION: Prox operator of the iso norm 
%   INPUT:  
%           t                = stepsize
%           (y2,y3)          = input vector
%           gamma            = de-noising parameter
%   OUTPUT:  gamma*prox( iso(y2,y3) )
%

Y23 = [y2, y3]; 
norm = [sqrt(y2.^2+y3.^2)]; %component-wise: sqrt(y2_k^2 + y3_k^2)
logic = norm>t; %vector of 1s and 0s to represent the constraint of the prox

%if norm > t: (y2_k y3_k) - t*(y2_k y3_k)/norm, aka divide each row [y2_k y3_k] by its norm

%if otherwise: 0 (due to logic vector)

%multiply in gamma bc this is prox of g(y)
output = gamma*(logic.*(Y23 - t*Y23./norm) );


end