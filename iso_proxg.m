function [prox2,prox3] = iso_proxg(y2,y3, t, gamma)
%DESCRIPTION: Prox operator of the iso norm 
%   INPUT:           
%           y2,y3            = y(:,:,2) and y(:,:,3)
%           gamma            = de-noising parameter (step-size)
%   OUTPUT:  prox( iso(y2,y3) )

 
    alpha = 1-gamma./max( (y2.^2+y3.^2).^(1/2) , gamma ); 

    prox2 = y2.*alpha;
    prox3 = y3.*alpha;

end