%% Initialization
clc;
clear;
close all;

%% Set basic parameters
[bits_original,x,y] = fImageSource(); % Read the original image

load X_task_B
X= X_task_B;
L_X = length(X);

load GoldSeq
GS = GoldSeq1;
L_GS = length(GS);

delay = [0 5 11 18 24 29];

x_1_delay = [X(delay(1)+1:end) zeros(1,delay(1))];
x_2_delay = [X(delay(2)+1:end) zeros(1,delay(2))];
x_3_delay = [X(delay(3)+1:end) zeros(1,delay(3))];
x_4_delay = [X(delay(4)+1:end) zeros(1,delay(4))];
x_5_delay = [X(delay(5)+1:end) zeros(1,delay(5))];
x_6_delay = [X(delay(6)+1:end) zeros(1,delay(6))];

x_1_r = GS*reshape(x_1_delay,L_GS,[]);
x_2_r = GS*reshape(x_2_delay,L_GS,[]);
x_3_r = GS*reshape(x_3_delay,L_GS,[]);
x_4_r = GS*reshape(x_4_delay,L_GS,[]);
x_5_r = GS*reshape(x_5_delay,L_GS,[]);
x_6_r = GS*reshape(x_6_delay,L_GS,[]);

beta = [0.2 0.5 0.9 0.5 0.4 0.6];
x_rake = beta*[x_1_r; x_2_r; x_3_r; x_4_r; x_5_r; x_6_r];
x_no = beta(1).*x_1_r;

%% Main Loop
L = length(bits_original);
BER_rake = 0;
BER_no = 0;

X_r = floor((sign(real(x_rake(1:L)))+1)/2);
X_n = floor((sign(real(x_no(1:L)))+1)/2);
for p = 1:L
    if X_r(1,p) ~= bits_original(p)
        BER_rake = BER_rake + 1;
    end
    if X_n(1,p) ~= bits_original(p)
        BER_no = BER_no + 1;
    end
end

%% Calculate BER
BER_rake = BER_rake/L
BER_no = BER_no/L

%% Regenerate the image
subplot(2,1,1)
[image_regenerate] = fImageSink(X_r,x,y);
f1 = imshow(image_regenerate);
title('RAKE')

subplot(2,1,2)
[image_regenerate] = fImageSink(X_n,x,y);
f2 = imshow(image_regenerate);
title('No RAKE')







