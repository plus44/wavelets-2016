function [SR_image] = Main_SR(thresh, s_or_h)
% *************************************************************************
% Wavelets and Applications Course - Dr. P.L. Dragotti
% MATLAB mini-project 'Sampling Signals with Finite Rate of Innovation'
% Exercice 6
% *************************************************************************
% 
% FOR STUDENTS
%
% This function runs a complete super-resolution algorithm based on the
% continuous moments analysis. The set of images is 'LR_Tiger_xx.tif' where
% xx is a number between 01 and 40. Each image is 64 x 64 pixels. The
% output super-resolved image 'SR_image' is 512 x 512 pixels.
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
   disp('Incorrect number of arguments to Main_SR.'); 
end

% Register images
[Tx_RGB, Ty_RGB]= ImageRegistration(thresh, s_or_h);

% Compute super-resolved image
[SR_image]= ImageFusion(Tx_RGB, Ty_RGB);