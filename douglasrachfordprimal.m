function [x] = douglasrachfordprimal(b, t, rho, init_vectors, problem, i, applyA_functions)
%DESCRIPTION: Primal Douglas-Rachford Splitting Algorithm
%   INPUT:  b                = blurred image
%           t                = stepsize
%           rho              = relaxation parameter
%           init_vectors     = the initial vectors (z_1, z_2)
%           problem          = specify the norm 
%   OUTPUT: x                = "deblurred" image 
%    

    applyA = applyA_functions{1};
    applyAT = applyA_functions{2};
    invertMatrix = applyA_functions{3};
    z_1 = init_vectors{1};
    z_2 = init_vectors{2};
    for k = 1:i.maxiter
        x = proxf(t,z_1); % proxf(z_1^k-1)
        if strcmp(problem, 'l1') == 1
           gamma = i.gammal1;
        else
           gamma = i.gammal2;
        end
        y = proxg(problem, b, z_2, t, gamma); % proxg(z_2^k-1)
        u = invertMatrix(2*x - z_1 + applyAT(2*y-z_2));
        v = applyA(u);
        z_1 = z_1 + rho*(u-x);
        z_2 = z_2 + rho*(v-y);
    end
    x = proxf(t,z_1);

end
