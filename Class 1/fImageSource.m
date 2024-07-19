function [bitsOut,x,y]=fImageSource(~)
filename = ('Flamingos.jpg');
i = importdata(filename);
images = imresize(i,0.3);
[x,y,z] = size(images);
bitsOut = dec2bin(images(:))';
bitsOut = double(bitsOut(:))-48;
