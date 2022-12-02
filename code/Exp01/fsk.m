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
A0=10;
A1=10;
fc1=5e3;
fc2=2*fc1;
cbits = 1.-bits;
y = bits.*(A1*sin(2*pi*fc1*t)) + cbits.*(A0*sin(2*pi*fc2*t));
% figure;
subplot(212);
plot(t, y);