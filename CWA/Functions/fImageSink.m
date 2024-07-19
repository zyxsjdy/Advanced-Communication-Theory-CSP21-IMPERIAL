% TimZ, MSc, 2021, Imperial College.
% 20/12/2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display the received image by converting bits back into R, B and G
% matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% photo_bit (Px1 Integers) = P demodulated bits of 1's and 0's
% b (Integer) = Number of bits in the image
% w (Integer) = Number of pixels in image in x dimension
% h (Integer) = Number of pixels in image in y dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% None
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fImageSink(photo_bit, b, w, h)

[~,n_signal] = size(photo_bit);  % number of signals

% Display the photo
figure;
for i = n_signal:-1:1
    bit_stream = photo_bit(1:b(i),i);  % Get rid of the added zeros
    dec_stream = bit2int(bit_stream,8);  % Convert it from binary stream to decimal stream
    photo = reshape(dec_stream,w(i),h(i),3);  % Reshape the decimal stream to the RGB format photo (x*y*3)
    subplot(n_signal,1,i);
    imshow(uint8(photo));  % Display the photo, "uint8(.)" means that we use a 8-bit number / 256 colours to represent the photo
end
end