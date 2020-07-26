% caromykle - August 2020
% Standard image import and formatting in MATLAB
clc;

A = imread('your_image.png');
A = im2double(A);

% Extracting one layer from the image
A = A(:,:,1);
