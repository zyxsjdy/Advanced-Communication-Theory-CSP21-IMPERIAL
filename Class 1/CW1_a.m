%% Initialization
clc;
clear;
close all;

%% Set basic parameters
load X_task_A
L = max(size(X_task_A));
for i = 1:6
    P_d(i) = sum(abs(X_task_A(i,:)).^2);
end
[~,I] = max((P_d-0.6*L)/(0.6*L));

beta = [0.2 0.5 0.9 0.5 0.4 0.6];
X_MRC = zeros(1,L);
X_SC = zeros(1,L);
[bits_original,x,y] = fImageSource(); % Read the original image
BER_MRC = 0;
BER_SC = 0;

%% Main Loop
for p = 1:L
    sd_1 = beta*X_task_A(:,p);
    sd_2 = X_task_A(:,p);
    
    X_MRC(1,p) = real(sd_1);
    if X_MRC(1,p) < 0
        X_MRC(1,p) = 1;
    else
        X_MRC(1,p) = 0;
    end
    if X_MRC(1,p) ~= bits_original(p)
        BER_MRC = BER_MRC + 1;
    end

    X_SC(1,p) = real(sd_2(I));
    if X_SC(1,p) < 0
        X_SC(1,p) = 1;
    else
        X_SC(1,p) = 0;
    end
    if X_SC(1,p) ~= bits_original(p)
        BER_SC = BER_SC + 1;
    end
end

%% Calculate BER
BER_MRC = BER_MRC/L
BER_SC = BER_SC/L

%% Regenerate the image
subplot(2,1,1)
bits_regenerate = X_MRC;
[image_regenerate] = fImageSink(bits_regenerate,x,y);
f1 = imshow(image_regenerate);
title('MRC')

subplot(2,1,2)
bits_regenerate = X_SC;
[image_regenerate] = fImageSink(bits_regenerate,x,y);
f2 = imshow(image_regenerate);
title('SC')

















