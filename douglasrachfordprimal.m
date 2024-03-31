function [x] = douglasrachfordprimal(b, t, rho, init_vectors, problem)
%DESCRIPTION: Primal Douglas-Rachford Splitting Algorithm
%   INPUT:  b                = blurred image
%           t                = stepsize
%           p                = relaxation parameter
%           init_vectors     = the initial vectors (z_1, z_2)
%           problem          = specify the norm 
%   OUTPUT: x                = "deblurred" image 
%
%    for k =1:i.maxiter
 %       x = prox_f(t,init_vectors[0]);
  %      if strcmp(problem,'l1') == 1
   %         gamma=i.gammal1;
    %    else
%            gamma = i.gammal2;
 %       end
  %      y = norm_prox(problem, b, init_vectors[1], t, rho, gamma);
   %     temp = 2*x - init_vectors[0] ...
    %end

end
