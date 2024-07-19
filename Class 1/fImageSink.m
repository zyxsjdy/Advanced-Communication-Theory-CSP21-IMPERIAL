function [image] = fImageSink(bitsIn,x,y)
type = 3; % 3 implies that this is a RGB image
Q = x*y*type;
img = reshape(bitsIn,length(bitsIn)/Q,Q);
image = zeros(1,Q);
for i = 1:Q
    image(:,i) = bin2dec(sprintf('%d',(img(:,i)')));
end
image = reshape(uint8(image),x,y,type);
