function [prox2,prox3] = iso_proxg(y2,y3, t, gamma)
%DESCRIPTION: Prox operator of the iso norm 
%   INPUT:  
%           t                = stepsize
%           y2,y3            = y(:,:,2) and y(:,:,3)
%           gamma            = de-noising parameter
%   OUTPUT:  gamma*prox( iso(y2,y3) )


    %if norm > t: (y2_k y3_k) - t*(y2_k y3_k)/norm, 
    % aka divide each row [y2_k y3_k] by its norm
    %if otherwise: 0 (due to logic vector)

    norm = (y2.^2+y3.^2).^(1/2); %component-wise: sqrt(y2_k^2 + y3_k^2)
    logic = norm>t; %bool vect represents the constraint of the prox+

    prox2 = gamma*(logic.*(y2 - t*y2./norm) ); 
    prox3 = gamma*(logic.*(y3 - t*y3./norm) );


end