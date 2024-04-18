%MAIN SCRIPT
%storing+blurring image:

    I = imread('True_Image.png'); %Store image
    figure('Name','image before deblurring') % Show initial image
    imshow(I,[])   
    I = double(I(:, :, 1));% Resize image (pixels between 0-1) 
    mn=min(I(:));
    I = I-mn;
    mx = max(I(:));
    I = I/mx;
    % can use this to resize image for faster computation
    %resizefactor = 0.1;
    %I = imresize(I, resizefactor);

    % Calculate the mean of the first column
    mean_first_column = mean(I(:, 1));
    % Calculate the mean of the first row
    mean_first_row = mean(I(1, :));
    
    m = (mean_first_row+mean_first_column) / 2;
    % Add padding around the image
    I = padarray(I,[7 7],m,'both');
    figure('Name','image before deblurring') % Show initial image
    imshow(I,[]) 

    % Generate blurred image
    noiseDensity = 0.5; 
    kernel = fspecial('gaussian', [15, 15], 5); 
    b = imfilter(I,kernel);
    %b = imnoise(b,'salt & pepper',noiseDensity);
    [numRows, numCols] = size(b);
    figure('Name','image after blurring') % Show blurred image
    imshow(b,[]) 

% default parameters:

    %common parameters
    i.maxiter = 500;
    i.gammal1 = 0;
    i.gammal2 = 0;
    %alg1
        % Set parameters for Alg1
        i.tprimaldr = 0.1;
        i.rhoprimaldr = 2.0;
        % Set initial vectors for Alg1
        z_1 = zeros(numRows,numCols);
        %z_1(1, 1)=1;
        z_2 = cat(3,z_1,z_1,z_1); % |z_2|=3n^2
        x_initAlg1 = {z_1, z_2};
    %alg2
        % Set parameters for Alg2
        i.tprimaldualdr = 2.0;
        i.rhoprimaldualdr = 1.049;
        % Set initial vectors for Alg2
        p = zeros(numRows, numCols);
        q = cat(3,p,p,p); % |q|=3n^2
        x_initAlg2 = {p,q};
    %alg 3
        % Set parameters for Alg3
        i.tadmm = 1.3;
        i.rhoadmm = 0.9;
        % Set initial vectors for Alg3
        u = zeros(numRows, numCols);
        y = cat(3,u,u,u); % |y|=3n^2
        w = zeros(numRows, numCols);
        z = cat(3,u,u,u); % |z|=3n^2
        x_initAlg3 = {u, y, w, z};

% Deblurring the image:

    x = optsolve('l1', 'admm', x_initAlg3, kernel, b, i);

    %figure('Name','image after deblurring') % Show deblurred image
    %imshow(x,[]) 
    
    pad_size = 7; % Example padding size, adjust as needed

    % Get the size of the original matrix
    [rows, cols] = size(x);
    
    % Define the range of rows and columns after removing padding
    row_range = (1 + pad_size):(rows - pad_size);
    col_range = (1 + pad_size):(cols - pad_size);
    cropped_matrix = x(row_range, col_range);
    figure('Name','image after deblurring') % Show deblurred image
    imshow(cropped_matrix,[]) 