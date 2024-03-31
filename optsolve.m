
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

    % ...

    % deblurring
    if strcmp(algorithm, 'douglasrachfordprimal') == 1
        x = douglasrachfordprimal(b, i.tprimaldr, i.rhoprimaldr, x_init, problem);

    elseif strcmp(algorithm ,'douglasrachfordprimaldual') == 1
        x = douglasrachfordprimal(b, i.tprimaldualdr, i.rhoprimaldualdr, x_init, problem);
        
    else 
        x = admm(b, i.tadmm, i.rhoadmm, x_init, problem);
    end 

end