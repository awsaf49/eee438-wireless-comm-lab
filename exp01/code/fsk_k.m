clc;
clear all;
close all;
N=5000;
A0=10;
A1=10;
Rb = 1000;
Tb = 1/Rb;
fc1=5e3;
fc2=2*fc1;
Fs = 20*fc1;
Ts = 1/Fs;
k = Fs*Tb;


% bit stream
bits = randi([0,1], 1, N);
bits = repelem(bits, k);
cbits = 1.-bits;
t= 0:Ts:N*Tb-Ts; %linspace(0, 0.5, N*k);

% modulation
kks = [1/8, 1/4, 1/2, 1, 2];
for kk = kks
    fc3 = fc1 + kk/Tb;
    y = bits.*(A1*sin(2*pi*fc1*t)) + cbits.*(A0*sin(2*pi*fc3*t));
    figure;
    [pxx f] = pwelch(y, 500, 200, 500, Fs);
    pxx = pxx/max(pxx);
    plot(f, pxx);
    xlim([0, 15000]);
    ylim([-0.1, 1.1]);
    title(sprintf("PSD @ k=%0.2f",kk))
    grid on;
end