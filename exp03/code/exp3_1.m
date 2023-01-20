clc;
clear all;
close all;

c = 3e8;
fc = 1.5e9;
d = 1:1:1000;
Gt=1;
Gr=1;
lam = c/fc;

%% free-space
num = Gt*Gr*(lam^2);
den = (4*pi*d).^2;
PL_fp = 10*log10(num./den);

semilogx(d, PL_fp);
grid on;
xlabel('log(distance)'), ylabel('PL free space (dB)');

%% long-distance
n = [2, 3, 6];
d0 = 100;

for i = 1:length(n)
    num = Gt*Gr*(lam^2);
    den = (4*pi*d0)^2;
    PL0_ld = -10*log10(num/den);
    PL_ld = PL0_ld + 10*n(i)*log10(d/d0);

    semilogx(d,PL_ld);
    grid on;
    xlabel('log(distance)'), ylabel('PL log-distance(dB)');
    hold on;
end

legend('n=2', 'n=3', 'n=6');

%% log-normal
num = Gt*Gr*(lam^2);
den = (4*pi*d0)^2;
sigma= 3;
n0=2;
% shadow=lognrand(0,3);
PL0_ln = 10*log10(num./den)+10*n0*log10(d/d0);
PL_ln= PL0_ln + sigma*randn(size(d));
% \lognrnd(0,3, size(d))
semilogx(d, PL_ln);
grid on;
xlabel('log(distance)'), ylabel('PL log normal shadow(dB)');

%% rayleigh
clc;
clear all;
close all;
N = 100000;
level = 50;
sigma = [1/sqrt(2), 1];
a= randn(1,N);
b= randn(1,N);
figure(1)

for sd=1:length(sigma)
    rl = sigma(sd)*(sqrt(a.^2+b.^2));
    [PDF, x] = hist(rl, level);
    mx = max(PDF);
    figure(1)
    plot(x, PDF)
    title('PDF');
    grid on; hold on;
    
    figure(2)
    plot(x, cumsum(PDF)/max(cumsum(PDF)))
    title('CDF');
    grid on; hold on;
    
end
figure(1)
legend('std=1/sqrt(2)', 'std=1');
figure(2)
legend('std=1/sqrt(2)', 'std=1');

%% rician

N = 100000;
level = 30;
sigma = 1/sqrt(2);
K_dB = [-40, 15];

figure(3);

for ki = 1:length(K_dB)
    K = 10^(K_dB(ki)/10);
    s = sqrt(K*2*sigma^2);

    e = s + sigma*randn(1, N);
    f = sigma*randn(1, N); 
    rc = sqrt(e.^2+f.^2);
    [PDF, x] = hist(rc, level);
    figure(3);
    plot(x, PDF);
    grid on;
    hold on;
    
    figure(4);
    plot(x, cumsum(PDF)/max(cumsum(PDF)))
    grid on;
    hold on;
end
figure(3);
legend('k = -40 dB', 'k = 15 dB');
figure(4);
legend('k = -40 dB', 'k = 15 dB');

%% Tx vs Rx Reyleigh
clc;
clear all;
close all;

N = 2000;
sigma = 1;
g=sigma*randn(1, N);
h=sigma*randn(1, N);
rl_sq=sqrt(g.^2+h.^2);
sq_sig = [ones(1, 400), -1*ones(1, 400), ones(1, 400),ones(1, 400),-ones(1, 400)];

rcv_sig = rl_sq.*sq_sig;

figure(5)
plot(rcv_sig);
hold on;
plot(sq_sig);
title('Reyleigh');

legend('Tx', 'Rx');

%% Tx vs Rx Rician

N = 2000;

K_dB = 15;
K = 10^(K_dB/10);
sigma = 1/sqrt(2);

mean1 = sqrt(K*2*sigma^2);

m = mean1 + sigma*randn(1, N); 
n = sigma*randn(1, N); 
sq_sig = [ones(1, 400), -1*ones(1, 400), ones(1, 400),ones(1, 400),-ones(1, 400)];
rcv_sig = sqrt(m.^2+n.^2).*sq_sig;
figure(6)
plot(rcv_sig);
hold on;
plot(sq_sig);
legend('Tx', 'Rx');
title('Rician')
