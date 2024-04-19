function [x] = admm(b, t, rho, init_vectors, problem, applyA_functions, i, eigArrys)
%DESCRIPTION: Alternating Direction Method of Multipliers Algorithm
%   INPUT:  b                = blurred image
%           t                = stepsize
%           p                = relaxation parameter
%           init_vectors     = the initial vectors (u_0, y_0, w_0, z_0)
%           problem          = specify the norm 
%   OUTPUT: x                = "deblurred" image 
%

    
    %setting up functions and variables
    [u_0, y_0, w_0, z_0] = init_vectors{:};
    applyA = applyA_functions{1};
    applyAT = applyA_functions{2};
    eigArry_K = eigArrys{1};
    eigArry_D1 = eigArrys{2};
    eigArry_D2 = eigArrys{3};
    eigArry_KTrans = eigArrys{4};
    eigArry_D1Trans = eigArrys{5};
    eigArry_D2Trans = eigArrys{6};
    % matrix and the eigenvalues of I + t*t*K^TK + t*t*D^TD; here t is the stepsizes
    eigValsMatT = ones(size(eigArry_K)) + eigArry_KTrans.*eigArry_K + eigArry_D1Trans.*eigArry_D1 + eigArry_D2Trans.*eigArry_D2;
    invertMatrixT = @(x) real(ifft2(fft2(x)./eigValsMatT)); 

    if strcmp(problem, 'l1') == 1
           gamma = i.gammal1;
    else
           gamma = i.gammal2;
    end

    for k=1:i.maxiter
        x_1 = invertMatrixT(u_0 + applyAT(y_0) - 1/t*( w_0 + applyAT(z_0) ) );
        u_1 = proxf(1/t, rho*x_1 + (1-rho)*u_0 + w_0/t  );
        y_1 = proxg(problem,b, rho*applyA(x_1) + (1-rho)*y_0 + z_0/t,1/t,gamma);
        w_1 = w_0 + t*(x_1-u_1);
        z_1 = z_0 + t*(applyA(x_1)-y_1);

        u_0 = u_1;
        y_0 = y_1;
        w_0 = w_1;
        z_0 = z_1;
    end

    x = invertMatrixT(u_1+ applyAT(y_1)- 1/t*(w_1+applyAT(z_1)) );
    

end
