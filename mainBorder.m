%MAIN SCRIPT
%storing+blurring image:
tic
    I = imread('mcgill.jpg'); %Store image
%     figure('Name','image before deblurring') % Show initial image
    %imshow(I,[])  
    I = rgb2gray(I);
    I = double(I(:, :, 1));% Resize image (pixels between 0-1) 
    mn=min(I(:));
    I = I-mn;
    mx = max(I(:));
    I = I/mx;
    % can use this to resize image for faster computation
    %resizefactor = 0.1;
    %I = imresize(I, resizefactor);

    I = padarray(I, [7 7], "symmetric");
%     figure('Name','image before deblurring') % Show initial image
    %imshow(I,[]) 

    % Generate blurred image
    noiseDensity = 0.5; 
    kernel = fspecial('gaussian', [15, 15], 5); 
    b = imfilter(I,kernel);
    b = imnoise(b,'salt & pepper',noiseDensity);
    [numRows, numCols] = size(b);
    figure('Name','image after blurring') % Show blurred image
    imshow(b,[]) 

% default parameters:

    %common parameters
    i.maxiter = 1;
%     i.gammal1 = 0.0076;
    i.gammal1 = 0.003;
    i.gammal2 = 0.0;
    %alg1
        % Set parameters for Alg1
        i.tprimaldr = 0.01;
        i.rhoprimaldr = 1.95;
        % Set initial vectors for Alg1
        z_1 = zeros(numRows,numCols);
        %z_1(1, 1)=1;
        z_2 = cat(3,z_1,z_1,z_1); % |z_2|=3n^2
        x_initAlg1 = {z_1, z_2};
    %alg2
        % Set parameters for Alg2
        i.tprimaldualdr = 10;
        i.rhoprimaldualdr = 1;
        % Set initial vectors for Alg2
        p = zeros(numRows, numCols);
        q = cat(3,p,p,p); % |q|=3n^2
        x_initAlg2 = {p,q};
    %alg 3
        % Set parameters for Alg3
        i.tadmm = 15;
        i.rhoadmm = 1.5;
        % Set initial vectors for Alg3
        u = zeros(numRows, numCols);
        y = cat(3,u,u,u); % |y|=3n^2
        w = zeros(numRows, numCols);
        z = cat(3,u,u,u); % |z|=3n^2
        x_initAlg3 = {u, y, w, z};

% Deblurring the image:

    x = optsolve('l1', 'douglasrachfordprimaldual', x_initAlg2, kernel, b, i);
    %figure('Name','image after deblurring') % Show deblurred image
    %imshow(x,[]) 
    
    pad_size = 7; % Example padding size, adjust as needed

    % Get the size of the original matrix
    [rows, cols] = size(x);
    
    % Define the range of rows and columns after removing padding
    row_range = (1 + pad_size):(rows - pad_size);
    col_range = (1 + pad_size):(cols - pad_size);
    cropped_matrix = x(row_range, col_range);
    toc
    L2SquaredError = norm(cropped_matrix - I)^2
    figure('Name','image after deblurring') % Show deblurred image
    imshow(cropped_matrix,[]) 