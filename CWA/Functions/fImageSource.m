% TimZ, MSc, 2021, Imperial College.
% 20/12/2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reads an image file with AxB pixels and produces a column vector of bits
% of length Q=AxBx3x8 where 3 represents the R, G and B matrices used to
% represent the image and 8 represents an 8 bit integer. If P>Q then
% the vector is padded at the bottom with zeros
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% photo_name (String) = The file name of the image
% b (Integer) = Number of bits to produce at the output - Should be greater
% than or equal to Q=AxBx3x8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% photo_bit (Px1 Integers) = P bits (1's and 0's) representing the image
% w (Integer) = Number of pixels in image in x dimension
% h (Integer) = Number of pixels in image in y dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [photo_bit, w, h] = fImageSource(photo_name, b)

% Initialize the output variable "photo_bit"
photo_bit = zeros(b,1);

% Read the photo and convert it to a bit stream
photo = imread(photo_name);
[w,h,~] = size(photo);  % Obtain the width and height of the photo
dec_stream = reshape(photo,[],1); % Reshape the photo matrix, change it to a decimal stream
bit_stream = int2bit(dec_stream,8); % Convert the decimal stream to a bit stream
photo_bit(1:length(bit_stream)) = bit_stream;
end