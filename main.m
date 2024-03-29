function [output] = prox_f(t, x)
%DESCRIPTION: Prox operator of f
%   INPUT:  t                = stepsize
%           x                = input vector
%   OUTPUT:  
%
end

function [output] = norm_prox(norm, b,x, t, rho, gamma)
%DESCRIPTION: Prox operator for l1 and l2 norms 
%   INPUT:  norm             = problem (either l1 or l2)
%           b                = blurred image
%           t                = stepsize
%           rho              = relaxation parameter
%           x                = input vector
%           gamma            = de-noising parameter
%   OUTPUT:  
%
end

function [output] = iso_prox(b,x, t, rho, gamma)
%DESCRIPTION: Prox operator of the iso norm 
%   INPUT:  b                = blurred image
%           t                = stepsize
%           rho              = relaxation parameter
%           x                = input vector
%           gamma            = de-noising parameter
%   OUTPUT:  
%
end

function [x] = douglasrachfordprimal(b, t, rho, init_vectors, problem)
%DESCRIPTION: Primal Douglas-Rachford Splitting Algorithm
%   INPUT:  b                = blurred image
%           t                = stepsize
%           p                = relaxation parameter
%           init_vectors     = the initial vectors (z_1, z_2)
%           problem          = specify the norm 
%   OUTPUT: x                = "deblurred" image 
%

end

function [x] = douglasrachfordprimaldual(b, t, rho, init_vectors, problem)
%DESCRIPTION: Primal-Dual Douglas-Rachford Splitting Algorithm
%   INPUT:  b                = blurred image
%           t                = stepsize
%           p                = relaxation parameter
%           init_vectors     = the initial vectors (p_0, q_0)
%           problem          = specify the norm 
%   OUTPUT: x                = "deblurred" image 
%

end

function [x] = admm(b, t, rho, init_vectors, problem)
%DESCRIPTION: Alternating Direction Method of Multipliers Algorithm
%   INPUT:  b                = blurred image
%           t                = stepsize
%           p                = relaxation parameter
%           init_vectors     = the initial vectors (x_0, u_0, y_0, w_0, z_0)
%           problem          = specify the norm 
%   OUTPUT: x                = "deblurred" image 
%

end


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


%Store image
I= imread('True_Image.png');

% Show initial image
figure('Name','image before deblurring')
imshow(I,[])

% Resize image by making pixels between 0 and 1 
I = double(I(:, :, 1));
mn=min(I(:));
I=I-mn;
mx = max(I(:));
I = I/mx;

% can use this to resize image for faster computation
%resizefactor = 0.1;
%I = imresize(I, resizefactor);

% Generate blurred image
noiseDensity = 0.5; 
kernel = fspecial('gaussian', [15, 15], 5); 
b = imfilter(I,kernel);
b = imnoise(b,'salt & pepper',noiseDensity);
[numRows, numCols] = size(b);

% Show blurred image
figure('Name','image after blurring')
imshow(b,[]) 

% Initialize commun parameters
i.maxiter = 500;
i.gammal1 = 0.049;
i.gammal2 = 0.049;

% Set parameters for Alg1
i.tprimaldr = 2.0;
i.rhoprimaldr = 0.1;
% Set initial vectors for Alg1
z_1 = zeros(numRows*numCols, 1);
z_2 = zeros(numRows*numCols, 1);
x_initAlg1 = [z_1, z_2];

% Set parameters for Alg2
i.tprimaldualdr = 2.0;
i.rhoprimaldualdr = 1.049;
% Set initial vectors for Alg2
p= zeros(numRows*numCols, 1);
q = zeros(numRows*numCols, 1);
x_initAlg2 = [p,q];

% Set parameters for Alg3
i.tadmm = 2.0;
i.rhoadmm = 1.049;
% Set initial vectors for Alg3
x= zeros(numRows*numCols, 1);
u = zeros(numRows*numCols, 1);
y = zeros(3*numRows*numCols, 1); % |y|=3n^2
w = zeros(numRows*numCols, 1);
z= zeros(3*numRows*numCols, 1); % |z|=3n^2
x_initAlg3 = {x, u, y, w, z};

% Deblurring the image
%x= optsolve('l1', 'douglasrachfordprimal', x_initAlg1, kernel, b, i);