clc;
clear all;
close all;

c = 3e8;
fc = 1.5e9;
d = 1:1:1000;
Gt=1;
Gr=1;
lam = c/fc;
d0 = 100;

%% free-space
num = Gt*Gr*(lam^2);
den = (4*pi*d).^2;
PL_fp = -10*log10(num./den);

semilogx(d, PL_fp);
grid on;
xlabel('log(distance)'), ylabel('Path Loss (dB)');
title('Free Space Model')

%% long-distance
n = [2, 3, 6];

for i = 1:length(n)
    num = Gt*Gr*(lam^2);
    den = (4*pi*d0)^2;
    PL0_ld = -10*log10(num/den);
    PL_ld = PL0_ld + 10*n(i)*log10(d/d0);

    semilogx(d,PL_ld);
    grid on;
    xlabel('log(distance)'), ylabel('Path Loss (dB)');
    title('Long-distance model')
    hold on;
end

legend('n=2', 'n=3', 'n=6');

%% log-normal shadowing
num = Gt*Gr*(lam^2);
den = (4*pi*d0)^2;
sigma= 3;
n0=2;
PL0_ln = 10*log10(num./den)+10*n0*log10(d/d0);
PL_ln= PL0_ln + sigma*randn(size(d));
semilogx(d, PL_ln);
grid on;
xlabel('log(distance)'), ylabel('Path Loss (dB)');
title('Log-normal shadowing')

%% rayleigh
clc;
clear all;
close all;
N = 100000;
level = 50;
sigma = [1/sqrt(2), 1, 1.5];
a= randn(1,N);
b= randn(1,N);
figure(1)

for sd=1:length(sigma)
    rl = sigma(sd)*(sqrt(a.^2+b.^2));
    [PDF, x] = hist(rl, level);
    mx = max(PDF);
    subplot(211)
    plot(x, PDF)
    title('Rayleigh Distribution (PDF)');
    xlabel('x'); ylabel('Occurance')
    grid on; hold on;
    
    subplot(212)
    plot(x, cumsum(PDF)/max(cumsum(PDF)))
    title('Rayleigh Distribution (CDF)');
    xlabel('x'); ylabel('Normalized Occurance')
    grid on; hold on;
    
end
subplot(211)
legend('std=1/sqrt(2)', 'std=1', 'std=1.5');
subplot(212)
legend('std=1/sqrt(2)', 'std=1', 'std=1.5');

%% rician

N = 100000;
level = 30;
sigma = 1/sqrt(2);
K_dB = [-40, 15, 10];

figure(3);

for ki = 1:length(K_dB)
    K = 10^(K_dB(ki)/10);
    s = sqrt(K*2*sigma^2);

    e = s + sigma*randn(1, N);
    f = sigma*randn(1, N); 
    rc = sqrt(e.^2+f.^2);
    [PDF, x] = hist(rc, level);
    subplot(211)
    plot(x, PDF);
    title('Rician Distribution (PDF)');
    xlabel('x'); ylabel('Occurance')
    grid on;
    hold on;
    
    subplot(212)
    plot(x, cumsum(PDF)/max(cumsum(PDF)))
    title('Rician Distribution (CDF)');
    xlabel('x'); ylabel('Normalizeed Occurance')
    grid on;
    hold on;
end
subplot(211)
legend('k = -40 dB', 'k = 15 dB', 'k = 5 dB');
subplot(212)
legend('k = -40 dB', 'k = 15 dB', 'k = 5 dB');

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

subplot(211)
plot(rcv_sig, 'LineWidth',1.3);
hold on; grid on;
plot(sq_sig, 'LineWidth',1.4);
hold off;
title('Reyleigh Distribution');
xlabel('x'); ylabel('Amp.')
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

subplot(212)
plot(rcv_sig, 'LineWidth',1.3);
hold on; grid on;
plot(sq_sig, 'LineWidth',1.4);
hold off; grid on;
title('Rician Distribution');
xlabel('x'); ylabel('Amp.')
legend('Tx', 'Rx');
