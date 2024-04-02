function [y] = proxg(norm, b,y, t, gamma)
%DESCRIPTION: Prox operator on g(y). first we find proxl1/l2 for y1 then we
%find iso prox for y2,y3. outputting the 3n x n matrix y(:,:,:) 
%   INPUT:  norm             = problem (either l1 or l2)
%           b                = blurred image for l1/l2
%           t                = stepsize
%           y                = input matrix 3nxn (concatenated y1,y2,y3)
%           gamma            = de-noising parameter for iso
%   OUTPUT: y = prox(t,g)

  
    %finding l1 or l2 prox
    if strcmp(norm,'l1') ==1
        %insert l1 prox g
        y(:,:,1) = proxl1(t,b,y(:,:,1));
    else
        %insert l2 prox g
        y(:,:,1) = l2prox(y(:,:,1),b,t);
    end

    %find iso prox
    [y(:,:,2), y(:,:,3)] = iso_proxg(y(:,:,2),y(:,:,3),t,gamma);

end
