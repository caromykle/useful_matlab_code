% caromykle - August 2020
% Standard image import and formatting in MATLAB
clc;

A = imread('your_image.png');
A = im2double(A);

% Extracting one layer from the image
A = A(:,:,1);

% Retrieving the size of the image stored in variables m and n
[m,n] = size(A);

% OR reformat image to arbitrary size [p,q] as follows
p = 200;
q = 200;
A = imresize(A,[p,q]);

% More compact format:
B = imread('your_image.png');
B = im2double(imresize(B(:,:,1),[p,q]));



