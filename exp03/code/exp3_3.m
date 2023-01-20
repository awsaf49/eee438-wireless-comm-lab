clc;
clear all;
close all;
d = 1:1:125000;
bits = randi([0,1], 1, 500);
bits = repelem(bits, 250);
t=linspace(0, 0.5, 500*250);
figure(1);
subplot(211);
stairs(t,bits)
% xaxis([0, 0.005]);
A0=10;
A1=10;
fc=5e3;
ph1 = 0;
ph2 = pi;
cbits = 1.-bits;
y = bits.*(A1*sin(2*pi*fc*t+ph1)) + cbits.*(A0*sin(2*pi*fc*t+ph2));
% figure;
subplot(212);
plot(t, y);



sigma= 3;
n0=2;
% shadow=lognrand(0,3);
PL0_ln = 0;
PL_ln= PL0_ln + sigma*randn(size(d));
% PL_ln = abs(Ray_model(125000));
% \lognrnd(0,3, size(d))
% figure(2)
% semilogx(d, PL_ln);

figure(3)
subplot(211)
sig_mul_fad = y.*PL_ln;
plot(d, sig_mul_fad);

subplot(212)
sig_add_noise = awgn(sig_mul_fad,-40);
plot(d, sig_add_noise);

