%% Initialization
clc;
clear;
close all;

%% figure(1)
load Ex03_Signal_Samples.mat
nt = x;
L = length(nt);

est_mean = sum(nt)/L
est_var = (sum((real(nt)).^2 + (imag(nt)).^2))/2 / L

%% figure(2)
mag = abs(nt);
count_mag = zeros(1,20);
step_mag = max(mag)/20;
for i = 1:L
    if mag(i) <= step_mag*1
        count_mag(1) = count_mag(1) + 1;
    elseif mag(i) <= step_mag*2
        count_mag(2) = count_mag(2) + 1;
    elseif mag(i) <= step_mag*3
        count_mag(3) = count_mag(3) + 1;
    elseif mag(i) <= step_mag*4
        count_mag(4) = count_mag(4) + 1;
    elseif mag(i) <= step_mag*5
        count_mag(5) = count_mag(5) + 1;
    elseif mag(i) <= step_mag*6
        count_mag(6) = count_mag(6) + 1;
    elseif mag(i) <= step_mag*7
        count_mag(7) = count_mag(7) + 1;
    elseif mag(i) <= step_mag*8
        count_mag(8) = count_mag(8) + 1;
    elseif mag(i) <= step_mag*9
        count_mag(9) = count_mag(9) + 1;
    elseif mag(i) <= step_mag*10
        count_mag(10) = count_mag(10) + 1;
    elseif mag(i) <= step_mag*11
        count_mag(11) = count_mag(11) + 1;
    elseif mag(i) <= step_mag*12
        count_mag(12) = count_mag(12) + 1;
    elseif mag(i) <= step_mag*13
        count_mag(13) = count_mag(13) + 1;
    elseif mag(i) <= step_mag*14
        count_mag(14) = count_mag(14) + 1;
    elseif mag(i) <= step_mag*15
        count_mag(15) = count_mag(15) + 1;
    elseif mag(i) <= step_mag*16
        count_mag(16) = count_mag(16) + 1;
    elseif mag(i) <= step_mag*17
        count_mag(17) = count_mag(17) + 1;
    elseif mag(i) <= step_mag*18
        count_mag(18) = count_mag(18) + 1;
    elseif mag(i) <= step_mag*19
        count_mag(19) = count_mag(19) + 1;
    elseif mag(i) <= step_mag*20
        count_mag(20) = count_mag(20) + 1;
    end
end

syms x_r
fx = x_r/(est_var)*exp(-(x_r^2)/(2*est_var));
integral = double(int(fx,x_r,0,step_mag*20));

x_ray = 0:step_mag/10:step_mag*20;
y_ray = x_ray/(est_var).*exp(-(x_ray.^2)/(2*est_var));

count_mag = count_mag ./ ((1000*step_mag)/(integral));

%% figure(3)
pha = angle(nt);
for i = 1:L
    if pha(i)<0
        pha(i) = 2*pi+pha(i);
    end
end
count_pha = zeros(1,20);
step_pha = 2*pi/20;
for i = 1:L
    if pha(i) <= step_pha*1
        count_pha(1) = count_pha(1) + 1;
    elseif pha(i) <= step_pha*2
        count_pha(2) = count_pha(2) + 1;
    elseif pha(i) <= step_pha*3
        count_pha(3) = count_pha(3) + 1;
    elseif pha(i) <= step_pha*4
        count_pha(4) = count_pha(4) + 1;
    elseif pha(i) <= step_pha*5
        count_pha(5) = count_pha(5) + 1;
    elseif pha(i) <= step_pha*6
        count_pha(6) = count_pha(6) + 1;
    elseif pha(i) <= step_pha*7
        count_pha(7) = count_pha(7) + 1;
    elseif pha(i) <= step_pha*8
        count_pha(8) = count_pha(8) + 1;
    elseif pha(i) <= step_pha*9
        count_pha(9) = count_pha(9) + 1;
    elseif pha(i) <= step_pha*10
        count_pha(10) = count_pha(10) + 1;
    elseif pha(i) <= step_pha*11
        count_pha(11) = count_pha(11) + 1;
    elseif pha(i) <= step_pha*12
        count_pha(12) = count_pha(12) + 1;
    elseif pha(i) <= step_pha*13
        count_pha(13) = count_pha(13) + 1;
    elseif pha(i) <= step_pha*14
        count_pha(14) = count_pha(14) + 1;
    elseif pha(i) <= step_pha*15
        count_pha(15) = count_pha(15) + 1;
    elseif pha(i) <= step_pha*16
        count_pha(16) = count_pha(16) + 1;
    elseif pha(i) <= step_pha*17
        count_pha(17) = count_pha(17) + 1;
    elseif pha(i) <= step_pha*18
        count_pha(18) = count_pha(18) + 1;
    elseif pha(i) <= step_pha*19
        count_pha(19) = count_pha(19) + 1;
    elseif pha(i) <= step_pha*20
        count_pha(20) = count_pha(20) + 1;
    end
end

x_uni = 0:step_pha:step_pha*20;
y_uni = ones(1,21);

count_pha = count_pha ./ ((1000*step_pha)/(2*pi));

%% Plot
figure(1)
plot(mag)

figure(2)
bar(step_mag*0.5:step_mag:step_mag*19.5,count_mag,1)
hold on
plot(x_ray,y_ray)

figure(3)
bar(step_pha*0.5:step_pha:step_pha*19.5,count_pha,1)
hold on
plot(x_uni,y_uni)


