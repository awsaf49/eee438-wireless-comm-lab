clear all;
close all;
clc;

Ts=50*1e-9;
fs=1/Ts;
sig_tao=10*1e-9; 
channel_number=10000; 
nbin=128;
P_max=ceil(10*sig_tao/Ts);

sig_0_sq=(1-exp(-Ts/sig_tao))/(1-exp(-(P_max+1)*Ts/sig_tao));

for p=0:P_max
    sig_p_sq(p+1)=sig_0_sq*exp(-p*Ts/sig_tao);
end
for m=1:length(sig_p_sq)
    h(:,m) = Ray_model(channel_number).*sqrt(sig_p_sq(m));
    avg_pow_h(m)= mean(h(:,m).*conj(h(:,m)));
    
end
H=fft(h,nbin);
H0 = fft(h(1,:),nbin);
H = mean(H, 2)';
H_all = fft(h,nbin);

figure(1)
stem([0:length(sig_p_sq)-1],sig_p_sq,'o'); 
hold on
stem([0:length(sig_p_sq)-1],avg_pow_h,'r.');
xlabel('p'), ylabel('Average Channel Power');
legend('Ideal','Simulation'); 
title('IEEE 802.11 channel')
grid on;

figure(2)
plot([-nbin/2+1:nbin/2]/nbin*fs/1e6,10*log10(H.*conj(H)),'b-');
xlabel('Frequency MHz'), ylabel('Channel power dB')
title('Average of Channels');
grid on;

figure(3)
plot([-nbin/2+1:nbin/2]/nbin*fs/1e6,10*log10(H0.*conj(H0)),'g-');
xlabel('Frequency MHz'), ylabel('Channel power dB')
title('Single Channel');
grid on;

figure(4)
plot([-nbin/2+1:nbin/2]/nbin*fs/1e6,10*log10(H_all.*conj(H_all)),'-');
xlabel('Frequency MHz'), ylabel('Channel power dB')
title('All Channels');
legend('1', '2', '3', '4', '5', '6')
grid on;



