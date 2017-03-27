function [Tx_RGB, Ty_RGB]= ImageRegistration(thresh, s_or_h)
% *************************************************************************
% Wavelets and Applications Course - Dr. P.L. Dragotti
% MATLAB mini-project 'Sampling Signals with Finite Rate of Innovation'
% Exercice 6
% *************************************************************************
% 
% FOR STUDENTS
%
% This function registers the set of 40 low-resolution images
% 'LR_Tiger_xx.tif' and returns the shifts for each image and each layer
% Red, Green and Blue. The shifts are calculated relatively to the first
% image 'LR_Tiger_01.tif'. Each low-resolution image is 64 x64 pixels.
%
%
% OUTPUT:   Tx_RGB: horizontal shifts, a 40x3 matrix
%           Ty_RGB: vertical shifts, a 40x3 matrix
%
% NOTE: _Tx_RGB(1,:) = Ty_RGB(1,:) = (0 0 0) by definition.
%       _Tx_RGB(20,2) is the horizontal shift of the Green layer of the
%       20th image relatively to the Green layer of the firs image.
%
%
% OUTLINE OF THE ALGORITHM:
%
% 1.The first step is to compute the continuous moments m_00, m_01 and m_10
% of each low-resolution image using the .mat file called:
% PolynomialReproduction_coef.mat. This file contains three matrices
% 'Coef_0_0', 'Coef_1_0' and 'Coef_0_1' used to calculate the continuous
% moments.
%
% 2.The second step consists in calculating the barycenters of the Red,
% Green and Blue layers of the low-resolution images.
%
% 3.By computing the difference between the barycenters of corresponding 
% layers between two images, the horizontal and vertical shifts can be 
% retrieved for each layer.
%
%
% Author:   Loic Baboulaz
% Date:     August 2006
%
% Imperial College London
% *************************************************************************

if nargin == 0,
    thresh = 0.35;
    s_or_h = 's'; 
elseif nargin == 1,
    s_or_h = 's';
elseif nargin > 2,
   disp('Incorrect number of arguments to ImageRegistration.'); 
end

% Load the coefficients for polynomial reproduction
load('PolynomialReproduction_coef.mat','Coef_0_0','Coef_1_0','Coef_0_1');

% -------- include your code here -----------
Tx_RGB = zeros(40, 3);
Ty_RGB = zeros(40, 3);
bc_ref = zeros(2, 3);
psf = fspecial('gaussian', 10, 1); 

% % % % % Perform a reference iteration on LR_Tiger_01 % % % % % 
filename = strcat('LR_Tiger_01.tif');
img = double(imread(filename)); % Load the image
img = edgetaper(img, psf); % Filter the edges
img = wthresh(img, s_or_h, thresh); % Threshold the image values

% Go through each color channel
for c = 1:3,
    % Find the moments of the current color channel
    m_00 = sum(sum(Coef_0_0.*img(:,:,c)));
    m_01 = sum(sum(Coef_0_1.*img(:,:,c)));
    m_10 = sum(sum(Coef_1_0.*img(:,:,c)));

    % Set the reference barycenter of the current color channel
    bc_ref(:, c) = [m_10/m_00; m_01/m_00];
    
    % Set the relative translation of the reference image as 0
    Tx_RGB(1, c) = 0;
    Ty_RGB(1, c) = 0;
end

% % % % % Iterate on all LR_Tiger_xx images % % % % % 
for n = 2:40,
    filename = strcat('LR_Tiger_', num2str(n, '%02d'));
    filename = strcat(filename, '.tif');
    
    img = double(imread(filename))./255;
    img = edgetaper(img, psf);
    img = wthresh(img, s_or_h, thresh);
    
    for c = 1:3,
        m_00 = sum(sum(Coef_0_0.*img(:,:,c)));
        m_01 = sum(sum(Coef_0_1.*img(:,:,c)));
        m_10 = sum(sum(Coef_1_0.*img(:,:,c)));
        
        % Find the barycenter of the current color channel of this img.
        bc = [m_10/m_00; m_01/m_00];
        
        % Set the relative translation by subtracting the reference bc.
        Tx_RGB(n, c) = bc(1) - bc_ref(1, c);
        Ty_RGB(n, c) = bc(2) - bc_ref(2, c);
    end
end



