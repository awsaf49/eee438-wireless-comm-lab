clc;
clear all;
close all;
M = 4;
i = (0:M-1);
phiM = i*2*pi/M;          
nb = log2(M);               

%% Constellation
A = 5;
ph = pi/4;
contX = A*cos(phiM+ph);
contY = A*sin(phiM+ph);

figure(1);
scatter(contX, contY, 'x');
title('QPSK constellation');
set(gca,'DataAspectRatio',[1 1 1])
grid on

%% Modulation

N0 = 200;
N = N0*nb;      
Rb = 1000;                  
k = 40;                     
Tb = 1/Rb;                  
fc = 2*Rb;                  
Fs = k*Rb;                  
Ts = 1/Fs;                  
t = 0:Ts:N*Tb-Ts;           

bits0 = round(rand(1,N));

symbols0 = bi2de(fliplr(reshape(bits0, [nb, N/nb])'))';
gray = bin2gray(symbols0,'psk',M);
phases = gray*2*pi/M; 

phase_up = length(t)/length(phases);
bit_up = length(t)/length(bits0);

bits = reshape(repmat(bits0, bit_up, 1), [1, length(t)]);
symbols = reshape(repmat(symbols0, phase_up, 1), [1, length(t)]);
phi = reshape(repmat(phases, phase_up, 1), [1, length(t)]);

mod = A*sin(2*pi*fc*t + phi + pi/4);
    
plot_len = 400;
figure(2);

subplot(311);
stairs(bits(1:plot_len));
ylim([-0.1, 1.1]);

subplot(312);
stairs(symbols(1:plot_len));
ylim([-0.5, 4.5]);

subplot(313);
plot(mod(1:plot_len));

%% PSD


figure;
[pxx f] = pwelch(mod, 500, 200, 500, Fs);
pxx = pxx/max(pxx);
plot(f, pxx);
% xlim([0, 10000]);
% ylim([-0.1, 1]);


%% PSD Window
pxx_list = [];
for idx=0:N-1
    win = mod(idx*k + 1:(idx+1)*k);
    win_pad = [zeros(1, idx*k), win, zeros(1, N*k - (idx+1)*k)];
    [pxx1 f] = pwelch(win_pad, 100, 50, 500, Fs);
    pxx1 = pxx1/max(pxx1);
    pxx_list = [pxx_list pxx1];
end
pxx_final = mean(pxx_list, 2);
pxx_final = pxx_final/max(pxx_final);

figure;
plot(f, pxx_final);