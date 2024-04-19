%MAIN SCRIPT
%storing+blurring image:
tic
     I= imread('True_Image.png'); %Store image
     figure('Name','image before deblurring') % Show initial image
     imshow(I,[])   
     I = double(I(:, :, 1));% Resize image (pixels between 0-1) 
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
     figure('Name','image after blurring') % Show blurred image
     imshow(b,[]) 
 
 % default parameters:
 
     %common parameters
     i.maxiter = 500;
     i.gammal1 = 0.003;%admm: 0.0076
     i.gammal2 = 0.0076;%0.089
     %alg1
         % Set parameters for Alg1
         i.tprimaldr = 0.01; % gamma = 0.003
         i.rhoprimaldr = 1.95;
         % Set initial vectors for Alg1
         z_1 = zeros(numRows, numCols);
         z_2 = cat(3,z_1,z_1,z_1);% |z_2|=3n^2
         x_initAlg1 = {z_1, z_2};
     %alg2
         % Set parameters for Alg2
         i.tprimaldualdr = 10.0;
         i.rhoprimaldualdr = 1.049;
         % Set initial vectors for Alg2
         p = zeros(numRows, numCols);
         q = cat(3,p,p,p); % |q|=3n^2
         x_initAlg2 = {p,q};
     %alg 3
         % Set parameters for Alg3
         i.tadmm = 1.03;
         i.rhoadmm = 0.85;%1.85;
         % Set initial vectors for Alg3
         u = zeros(numRows, numCols);
         y = cat(3,u,u,u); % |y|=3n^2
         w = zeros(numRows, numCols);
         z = cat(3,u,u,u); % |z|=3n^2
         x_initAlg3 = {u, y, w, z};
 
 % Deblurring the image:
 
     x = optsolve('l1', 'douglasrachfordprimal', x_initAlg1, kernel, b, i);
     toc
     figure('Name','image after deblurring') % Show blurred image
     imshow(x,[]) 

%%%%%%
%MAIN SCRIPT
%storing+blurring image:
tic
     I= imread('True_Image.png'); %Store image
     figure('Name','image before deblurring') % Show initial image
     imshow(I,[])   
     I = double(I(:, :, 1));% Resize image (pixels between 0-1) 
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
     figure('Name','image after blurring') % Show blurred image
     imshow(b,[]) 
 
 % default parameters:
 
     %common parameters
     i.maxiter = 500;
     i.gammal1 = 0.003;%admm: 0.0076
     i.gammal2 = 0.0076;%0.089
     %alg1
         % Set parameters for Alg1
         i.tprimaldr = 0.01; % gamma = 0.003
         i.rhoprimaldr = 1.95;
         % Set initial vectors for Alg1
         z_1 = zeros(numRows, numCols);
         z_2 = cat(3,z_1,z_1,z_1);% |z_2|=3n^2
         x_initAlg1 = {z_1, z_2};
     %alg2
         % Set parameters for Alg2
         i.tprimaldualdr = 10.0;
         i.rhoprimaldualdr = 1.049;
         % Set initial vectors for Alg2
         p = zeros(numRows, numCols);
         q = cat(3,p,p,p); % |q|=3n^2
         x_initAlg2 = {p,q};
     %alg 3
         % Set parameters for Alg3
         i.tadmm = 1.03;
         i.rhoadmm = 0.85;%1.85;
         % Set initial vectors for Alg3
         u = zeros(numRows, numCols);
         y = cat(3,u,u,u); % |y|=3n^2
         w = zeros(numRows, numCols);
         z = cat(3,u,u,u); % |z|=3n^2
         x_initAlg3 = {u, y, w, z};
 
 % Deblurring the image:
 
     x = optsolve('l1', 'douglasrachfordprimal', x_initAlg1, kernel, b, i);
     toc
     figure('Name','image after deblurring') % Show blurred image
     imshow(x,[]) 
