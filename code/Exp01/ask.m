clc;
clear all;
close all;
bits = randi([0,1], 1, 500);
bits = repelem(bits, 250);
t=linspace(0, 0.5, 500*250);
figure;
subplot(211);
stairs(t,bits)
% xaxis([0, 0.005]);
A0=5;
A1=10;
fc=5e3;
cbits = 1.-bits;
y = bits.*(A1*sin(2*pi*fc*t)) + cbits.*(A0*sin(2*pi*fc*t));
% figure;
subplot(212);
plot(t, y);