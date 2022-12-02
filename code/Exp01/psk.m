clc;
clear all;
close all;
N=500;
A0=10;
A1=10;
ph1 = 0;
ph2 = pi;
fc=5e3;
Rb = 1000;
Tb = 1/Rb;
Fs = 20*fc;
Ts = 1/Fs;
k = Fs*Tb;


bits = randi([0,1], 1, N);
bits = repelem(bits, k);
t= 0:Ts:N*Tb-Ts; %linspace(0, 0.5, N*k);
figure;
subplot(211);
stairs(t,bits)
% xaxis([0, 0.005]);
cbits = 1.-bits;
y = bits.*(A1*sin(2*pi*fc*t+ph1)) + cbits.*(A0*sin(2*pi*fc*t+ph2));
% figure;
subplot(212);
plot(t, y);

figure;
[pxx f] = pwelch(y, 500, 200, 500, Fs);
pxx = pxx/max(pxx);
plot(f, pxx);
xlim([0, 10000]);
ylim([-0.1, 1]);