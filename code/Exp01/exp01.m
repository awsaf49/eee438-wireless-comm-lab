clc;
clear all;
close all;
N=500;
A0=5;
A1=10;
fc=5e3;
Rb = 1000;
Tb = 1/Rb;
Fs = 40*fc;
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
y = bits.*(A1*sin(2*pi*fc*t)) + cbits.*(A0*sin(2*pi*fc*t));
% figure;
subplot(212);
plot(t, y);

figure;
[pxx f] = pwelch(y, 500, 200, 500, Fs);
pxx = pxx/max(pxx);
plot(f, pxx);
xlim([0, 10000]);
ylim([-0.1, 1]);