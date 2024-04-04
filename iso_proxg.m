function [prox2,prox3] = iso_proxg(y2,y3, t, gamma)
%DESCRIPTION: Prox operator of the iso norm 
%   INPUT:           
%           y2,y3            = y(:,:,2) and y(:,:,3)
%           gamma            = de-noising parameter (step-size)
%   OUTPUT:  prox( iso(y2,y3) )

    
    if gamma == 0
        %if gamma is 0 just output 0 (for debugging purposes)
        prox2 = zeros(size(y2));
        prox3 = zeros(size(y3));
    else

        %actual prox calculation
        alpha = 1-gamma./max( (y2.^2+y3.^2).^(1/2) , gamma ); 
    
        prox2 = y2.*alpha;
        prox3 = y3.*alpha;

    end