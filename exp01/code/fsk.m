clc;
clear all;
close all;
N=5000;
A0=10;
A1=10;
fc1=5e3;
fc2=2*fc1;
Rb = 1000;
Tb = 1/Rb;
Fs = 20*fc1;
Ts = 1/Fs;
k = Fs*Tb;


bits = randi([0,1], 1, N);
bits = repelem(bits, k);
t= 0:Ts:N*Tb-Ts; %linspace(0, 0.5, N*k);


figure;
subplot(211);
stairs(t,bits)
xlim([0, 0.012]);
ylim([-0.2, 1.2]);
title("Bit Stream")
grid on;


cbits = 1.-bits;
y = bits.*(A1*sin(2*pi*fc1*t)) + cbits.*(A0*sin(2*pi*fc2*t));
subplot(212);
plot(t, y);
xlim([0, 0.012]);
ylim([-12, 12]);
title("Modulation")
grid on;

figure;
[pxx f] = pwelch(y, 500, 200, 500, Fs);
pxx = pxx/max(pxx);
plot(f, pxx);
xlim([0, 40000]);
ylim([-0.1, 1.1]);
title("PSD")
grid on;