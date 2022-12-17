clc;
clear all;
close all;
N=5000; % number of points
A0=5; % symbol for 0
A1=10; % symbol for 1
fc=5e3; % carrier frequence
Rb = 1000; % bit rate (num_bits in 1s)
Tb = 1/Rb; % bit duration (time of 1 bit)
Fs = 40*fc; % sampling frequence
Ts = 1/Fs; % sampling time
k = Fs*Tb; % number of samples in 1 bit


bits = randi([0,1], 1, N);
bits = repelem(bits, k); % replicate bit to generate actual signal
t= 0:Ts:N*Tb-Ts; %linspace(0, 0.5, N*k);

figure;

% bitstream
subplot(211);
stairs(t,bits);
xlim([0, 0.012]);
ylim([-0.2, 1.2]);
title("Bit Stream")
grid on;

% modulated
cbits = 1.-bits;
y = bits.*(A1*sin(2*pi*fc*t)) + cbits.*(A0*sin(2*pi*fc*t));
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
xlim([0, 80000]);
ylim([-0.1, 1.1]);
title("PSD")
grid on;