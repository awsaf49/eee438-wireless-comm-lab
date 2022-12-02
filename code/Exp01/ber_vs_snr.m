clc;
clear all;
close all;
N=1e5;
A0=5;
A1=10;
ph1 = 0;
ph2 = pi;
fc1=5e3;
fc2=2*fc1;
Rb = 1000;
Tb = 1/Rb;
Fs = 20*fc1;
Ts = 1/Fs;
k = Fs*Tb;
vth_ask = (A0 + A1)/4;
vth_psk = 0;


b = randi([0,1], 1, N);
bits = repelem(b, k);
t= 0:Ts:N*Tb-Ts; %linspace(0, 0.5, N*k);
figure;
subplot(211);
stairs(t,bits)
% xaxis([0, 0.005]);
cbits = 1.-bits;
mc1 = sin(2*pi*fc1*t+ph1);
mc2 = sin(2*pi*fc1*t+ph2);
mc3 = sin(2*pi*fc2*t+ph1);
y_ask = bits.*(A1*mc1) + cbits.*(A0*mc1);
y_psk = bits.*(A1*mc2) + cbits.*(A0*mc1);
y_fsk = bits.*(A1*mc3) + cbits.*(A0*mc1);
% figure;
subplot(212);
plot(t, y_ask);

SNR = [-15, -10, -5, 0, 5];

yint_ask = [];
yint_psk = [];
yint_fsk = [];

yber_ask = [];
yber_psk = [];
yber_fsk = [];

count = 0;
for i=SNR
    count=count+1;
    yn_ask = awgn(y_ask, i, 'measured');
    yn_ask = yn_ask.*mc1;
    
    yn_psk = awgn(y_psk, i, 'measured');
    yn_psk = yn_psk.*mc1;
    
    yn_fsk = awgn(y_fsk, i, 'measured');
    yn_fsk1 = yn_fsk.*mc1;
    yn_fsk2 = yn_fsk.*mc3;
    
    
    for j=1:N
        yint_ask(j) = sum(yn_ask(k*(j-1)+1:k*j))/(k);
        yint_psk(j) = sum(yn_psk(k*(j-1)+1:k*j))/(k);
        yint_fsk1(j) = sum(yn_fsk1(k*(j-1)+1:k*j))/(k);
        yint_fsk2(j) = sum(yn_fsk2(k*(j-1)+1:k*j))/(k);
    end
    yd_ask = double(yint_ask>vth_ask);
    yd_psk = double(yint_psk<vth_psk);
    yd_fsk = double(yint_fsk1<yint_fsk2);
    
    %ber
    yber_ask(count) = sum(abs(yd_ask-b))/N;
    yber_psk(count) = sum(abs(yd_psk-b))/N;
    yber_fsk(count) = sum(abs(yd_fsk-b))/N;
end

figure;
semilogy(SNR, yber_ask,'r'); hold on;
semilogy(SNR, yber_psk,'b'); hold on;
semilogy(SNR, yber_fsk,'g'); hold on;
legend('bask','bpsk','bfsk');
% ylim([0,1]);