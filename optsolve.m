
function [x]= optsolve(problem,algorithm,x_init, kernel, b, i)
%DESCRIPTION: 
%   INPUT:  problem    = specify the norm use in the computation (either l1 or l2 norm) 
%           algorithm  = specify the algorithm that will be use for
%                        computation
%           x_init     = the initial vectors for the specified algorithm
%           kernel     = kernel used to blur the image
%           b          = blurred image
%           i          = stucture containing the input parameters
%   OUTPUT: x          = "deblurred" image 
%
    
    %Constructing the K and D matrices
    [numRows, numCols] = size(b); %numRows = m, numCols = n

    %computes the numRow x numCol matrix of the eigenvalues for K and D1 and
    %D2; Here D1 = I oplus D1 in the paper and D2 = D1 oplus I.
    eigArry_K = eigValsForPeriodicConvOp(kernel, numRows, numCols);
    eigArry_D1 = eigValsForPeriodicConvOp([-1,1], numRows, numCols);
    eigArry_D2 = eigValsForPeriodicConvOp([-1,1], numRows, numCols);

    %computes numRow x numCol matrix of the eigenvalues for K^T and D1^T and
    %D2^T;
    eigArry_KTrans = conj(eigArry_K);
    eigArry_D1Trans = conj(eigArry_D1);
    eigArry_D2Trans = conj(eigArry_D2);

    %Functions which compute Kx, D1x, D2x, Dxt, K^Tx, D1^Tx, D2^Tx, and D^Ty.
    %Note for all the x functions, the input x is in R^(m x n) and outputs into
    %R^(m x n) except for D which outputs into 2 concat. R^(m x n) matrices;
    %For D^Ty, y is two m x n matrices concatanated and outputs into R^(m x n)
    applyK = @(x) applyPeriodicConv2D(x, eigArry_K);
    applyD1 = @(x) applyPeriodicConv2D(x, eigArry_D1);
    applyD2 = @(x) applyPeriodicConv2D(x, eigArry_D2);
    
    applyKTrans = @(x) applyPeriodicConv2D(x, eigArry_KTrans);
    applyD1Trans = @(x) applyPeriodicConv2D(x, eigArry_D1Trans);
    applyD2Trans = @(x) applyPeriodicConv2D(x, eigArry_D2Trans);
    
    applyD = @(x) cat(3, applyD1(x), applyD2(x));
    
    applyDTrans = @(y) applyD1Trans(y(:,:,1)) + applyD2Trans(y(:, :, 2));

    %Functions which compute Ax and A^Ty
    %applyA takes 1 matrix and gives 3 concatenated matrices
    %applyAT takes 3 concatenated matrices and gives 1 matrix
    %Treat all (3n x n) matrices as (n x n x 3) for consistency and
        %simplicity
    %These functions can be used in algorithms
    applyA = @(x) cat(3, applyK(x), applyD(x));
    applyAT = @(y) applyKTrans(y(:,:,1)) + applyDTrans(y(:,:,2:3));
    % Function which computes the (I + K^TK + D^TD)x where x in R^(m x n)
    % matrix and the eigenvalues of I + t*t*K^TK + t*t*D^TD; here t is the
    % stepsizes
    applyMat = @(x) x + applyKTrans(applyK(x)) + applyDTrans(applyD(x));
    eigValsMat = ones(numRows, numCols) + eigArry_KTrans.*eigArry_K + eigArry_D1Trans.*eigArry_D1...
        + eigArry_D2Trans.*eigArry_D2;
    
    %R^(m x n) Computing (I + K^T*K + D^T*D)^(-1)*x
    invertMatrix = @(x) ifft2(fft2(x)./eigValsMat); 
    applyA_functions = {applyA, applyAT, invertMatrix};


    % deblurring
    if strcmp(algorithm, 'douglasrachfordprimal') == 1
        x = douglasrachfordprimal(b, i.tprimaldr, i.rhoprimaldr, x_init, problem, i, applyA_functions);

    elseif strcmp(algorithm ,'douglasrachfordprimaldual') == 1
        x = douglasrachfordprimal(b, i.tprimaldualdr, i.rhoprimaldualdr, x_init, problem, i);
        
    else 
        x = admm(b, i.tadmm, i.rhoadmm, x_init, problem, i);
    end 

end