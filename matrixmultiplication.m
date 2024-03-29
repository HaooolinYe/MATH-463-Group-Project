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

%Constructing the K and D matrices
[numRows, numCols] = size(b); %numRows = m, numCols = n

%computes the numRow x numCol matrix of the eigenvalues for K and D1 and
%D2; Here D1 = I oplus D1 in the paper and D2 = D1 oplus I.
eigArry_K = eigValsForPeriodicConvOp(kernel, numRows, numCols);
eigArry_D1 = eigValsForPeriodicConvOp([-1,1]', numRows, numCols);
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
t = i.tprimaldr; % Need to change this for the various algorithms you are applying
applyMat = @(x) x + applyKTrans(applyK(x)) + applyDTrans(applyD(x));
eigValsMat = ones(numRows, numCols) + t*t*eigArry_KTrans.*eigArry_K + t*t*eigArry_D1Trans.*eigArry_D1...
    + t*t*eigArry_D2Trans.*eigArry_D2;

%R^(m x n) Computing (I + K^T*K + D^T*D)^(-1)*x
invertMatrix = @(x) ifft2(fft2(x)./eigValsMat); 