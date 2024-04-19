function [x] = douglasrachfordprimaldual(b, init_vectors, problem, i, applyA_functions, eigArrys)
%DESCRIPTION: Primal-Dual Douglas-Rachford Splitting Algorithm
%   INPUT:  b                = blurred image
%           t                = stepsize
%           p                = relaxation parameter
%           init_vectors     = the initial vectors (p_0, q_0)
%           problem          = specify the norm 
%   OUTPUT: x                = "deblurred" image 
%
    %Initial parameters & necessary functions
    t = i.tprimaldualdr;
    rho = i.rhoprimaldualdr;
    if (strcmp(problem, 'l1')==1)
        gamma = i.gammal1;
    else
        gamma = i.gammal2;
    end
    p = init_vectors{1};
    q = init_vectors{2};
    applyA = applyA_functions{1};
    applyAT = applyA_functions{2};

    eigArry_K = eigArrys{1};
    eigArry_D1 = eigArrys{2};
    eigArry_D2 = eigArrys{3};
    eigArry_KTrans = eigArrys{4};
    eigArry_D1Trans = eigArrys{5};
    eigArry_D2Trans = eigArrys{6};

    % matrix and the eigenvalues of I + t*t*K^TK + t*t*D^TD; here t is the
    % stepsizes
    eigValsMatT = ones(size(eigArry_K)) + t*t*eigArry_KTrans.*eigArry_K + t*t*eigArry_D1Trans.*eigArry_D1 + t*t*eigArry_D2Trans.*eigArry_D2;
    invertMatrixT = @(x) real(ifft2(fft2(x)./eigValsMatT)); %function to apply (I+t^2A^TA)^-1
    for k = 1:i.maxiter

        %Resolvent of A
        x = proxf(t,p);     
        z = q - t*proxg(problem,b,q/t,1/t,gamma); %using moreau decomp. thm. to compute prox of g*

        %Resolvent of B (using computation from project description)
        matrixOne = 2*x - p;
        matrixTwo = 2*z - q;
        inverted = invertMatrixT(matrixOne-t*applyAT(matrixTwo));
        w = inverted;
        v = matrixTwo + t*applyA(inverted);

        p = p + rho*(w - x);
        q = q + rho*(v - z);
    end
    x = proxf(t, p);
end