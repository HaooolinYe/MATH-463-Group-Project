function [ eigValArry ] = eigValsForPeriodicConvOp( filter, numRows, numCols  )
%DESCRIPTION: Computes the eigenvalues of the 2D DFT of the convolution kernel;
%(Note this is an array of eigenvalues because numCols can be larger than).
%
%   INPUT:      filter  = correlation kernel (e.g. for Gaussian
%                       filter = fspecial('gaussian');
%               numRows = scalar; number of rows in the blurred image (b)
%               numCols = scalar; number of columns in the blurred image
%                       (b)
%   OUTPUT:     eigValArry = numRows x numCol matrix containing the
%                           eigenvalues of the convolution kernel

    % Constructing the impulse: customary to put this in the upper left hand
    % corner pixel
    a = zeros(numRows, numCols);
    a(1,1) = 1;

    %Impulse Response from the given kernel; 'circular' for periodic boundary
    %conditions
    Ra = imfilter(a, filter, 'circular');

    %Fourier transform of the impulse response (hence the hat)
    RaHat = fft2(Ra);

    eigValArry = RaHat;
end

function [ xtw ] = K_applyXTrans( w, x, kernelsize  )
%UNTITLED Summary of this function goes here
%   w \in \R^(256 x 256);
%   x \in \R^(256 x 256);
%   kernelsize = length of one side;

%output xtw = X^Tw

% (X^Tw)_(ij) = < w, conv2(x, impulse at the ij entry, 'same') >

    im = zeros(kernelsize);
    xtw = im;

    for i = 1:kernelsize
        for j = 1:kernelsize
            shift_im = im;
            shift_im(i,j) = 1;
            xtw(i,j) = sum(sum(w.*conv2(x, shift_im, 'same')));
        end
    end
    
end

function [ out ] = applyPeriodicConv2D( x, eigValArr )
%DESCRIPTION: For a given "unblurred" image x and the eigenvalue array for
%the blurring kernel, computes the "blurred image" (e.g. Kx and Dx in the
%paper)
%   INPUT:  x           = m x n matrix representing the image
%           eigValArry  = m x n representing the eigenvalues of the 2D DFT
%                       of the convolution kernel
%   OUTPUT: out         = m x n the "blurred" image (i.e. Kx)
%
%MATH: Observe that K = Q^H eigValArry Q where Q is essentially the
%discrete fourier transform and Q^H is the inverse fourier transform. This
%perform Kx which reduces to this ifft(eigValArry.*fft(x)).

    out = ifft2(eigValArr.*fft2(x));
end

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
%           gamma            = de-noising 
%   OUTPUT:  
%
end

function [output] = iso_prox(b,x, t, rho, gamma)
%DESCRIPTION: Prox operator of the iso norm 
%   INPUT:  b                = blurred image
%           t                = stepsize
%           rho              = relaxation parameter
%           x                = input vector
%           gamma            = de-noising 
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
    % matrix multiplication 
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

    % Function which computes the (I + K^TK + D^TD)x where x in R^(m x n)
    % matrix and the eigenvalues of I + t*t*K^TK + t*t*D^TD; here t is the
    % stepsizes
    if 'algorithm' == 'douglasrachfordprimal'
        t=i.tprimaldr
    elseif 'algorithm' == 'douglasrachfordprimaldual'
        t=i.tprimaldualdr
    else
        t=i.tadmm
    end
    %t = i.tprimaldr; % Need to change this for the various algorithms you are applying
    applyMat = @(x) x + applyKTrans(applyK(x)) + applyDTrans(applyD(x));
    eigValsMat = ones(numRows, numCols) + t*t*eigArry_KTrans.*eigArry_K + t*t*eigArry_D1Trans.*eigArry_D1...
        + t*t*eigArry_D2Trans.*eigArry_D2;

    %R^(m x n) Computing (I + K^T*K + D^T*D)^(-1)*x
    invertMatrix = @(x) ifft2(fft2(x)./eigValsMat); 
   
    %end matrix multiplication

    % deblurring
    if 'algorithm' == 'douglasrachfordprimal'
        x = douglasrachfordprimal(b, i.tprimaldr, i.rhoprimaldr, x_init, problem);

    elseif 'algorithm' == 'douglasrachfordprimaldual'
        x = douglasrachfordprimal(b, i.tprimaldualdr, i.rhoprimaldualdr, x_init, problem);
        
    else 
        x = admm(b, i.tadmm, i.rhoadmm, x_init, problem);
    end 

end


%Storing image in matrix 
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

% Resize image for faster computation
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
z_1 = zeros(numRows*numCols, 1);
z_2 = zeros(numRows*numCols, 1);
x_initAlg1 = [z_1, z_2];

% Set parameters for Alg2
i.tprimaldualdr = 2.0;
i.rhoprimaldualdr = 1.049;
p= zeros(numRows*numCols, 1);
q = zeros(numRows*numCols, 1);
x_initAlg2 = [p,q];

% Set parameters for Alg3
i.tadmm = 2.0;
i.rhoadmm = 1.049;
x= zeros(numRows*numCols, 1);
u = zeros(numRows*numCols, 1);
y = zeros(numRows*numCols, 1);
w = zeros(numRows*numCols, 1);
z= zeros(numRows*numCols, 1);
x_initAlg3 = [x, u, y, w, z];

% Deblur the image
%x= optsolve('problem', 'algorithm', x_init, kernel, b, i);