clc;
clear all;
close all;

%% task 1(i)
clc;
clear all;
close all;

Nt = 2;
Nr = 2;
N0 = 1;

s = 2*randi([0 1], Nt, 1) - 1;

noise = sqrt(N0/2) * (randn(Nr, 1) + randn(Nr, 1)); 
H = [0.4, 0.7; 1, 0.2];
y = H*s + noise;

figure(1);
plot(real(y), imag(y), 'o');
xlabel('Re axis');
ylabel('Im axis');

%% task 1(ii)

clc;
clear all;
close all;

Nt = 2;
Nr = 2;
N0 = 0.1; 

bits = [1, 1, -1, -1; 1, -1, 1, -1];

H = [0.4, 0.7; 1, 0.2];

s = 2*randi([0 1], Nt, 1) - 1;

noise = sqrt(N0/2) * (randn(Nr, 1) + randn(Nr, 1)); 

y = H*s + noise;
disp(['Transmitted Symbols: ', num2str(s')]);
disp(['Received Symbols: ', num2str(y')]);

mu = sum(abs(y - H*bits).^2);

[min_mu, idx] = min(mu);

detected = bits(:, idx);
disp(['Detected Symbols with MLD: ', num2str(detected')]);

%% task 1 (iii)

clc;
close all;
clear all;

Nt = 2;
Nr_vals = [2, 2];
bits = [1, 1, -1, -1; 1, -1, 1, -1];



sim_len = 1e4;
figure(1)

for sim_no  = 1:length(Nr_vals)
    Nr = Nr_vals(sim_no);
    SNR_values = -15:10;
    BER = zeros(size(SNR_values));
    H = [0.4, 0.7; 1, 0.2];
    for idx=1:length(SNR_values)
        s = 2*randi([0 1], Nt, sim_len) - 1;
        SNR = SNR_values(idx);
        N0 = 10^(-SNR/10);
        noise = sqrt(N0/2)*(randn(Nr, sim_len) + 1i*randn(Nr, sim_len));
        y = H*s + noise;
        y = repmat(y, 1, 1, 4);

        y1 = H*bits;
        y1 = repmat(y1, 1, 1, sim_len);
        y1 = permute(y1, [1, 3, 2]);

        mu = sum(abs(y - y1).^2);
        
        [min_mu, idx2] = min(mu, [], 3);
        
        detected = bits(:, idx2);

        se = (s ~= detected);
        ne = sum(se(:))/Nt;
        BER(idx) = ne/(Nt*sim_len);
    end
    plot(SNR_values, BER);
    xlabel('SNR'); ylabel('BER');
    
end

%% task2 (i)

clc;
clear all;
close all;


Nt_list = [1, 2];
Nr_list = [1, 2];

SNR_db = -15:0.01:10;
N0 = 10.^(-SNR_db/10);
N0_inv = 1./N0;
ch_name = [];

for Nt = Nt_list
    for Nr = Nr_list
        H = (randn(Nr, Nt) + 1i*randn(Nr, Nt))/sqrt(2);
        lam = eig(H*H');
        C = sum(log2(1 + lam*N0_inv/Nt));
        plot(SNR_db, C);
        hold on;
        ch_name = [ch_name, string(Nt)+'x'+string(Nr)];
    end
end
legend(ch_name);
xlabel('SNR');
ylabel('Capacity');

%% task2 (ii)
clc;
clear all;
close all;
SNR_db = -15:0.01:10;
N0 = 10.^(-SNR_db/10);
N0_inv = 1./N0;

H_2i2o = [1, 0.7; 0.4, 0.9];
H_siso=10;
H_2iso=[2, 0.7];
H_si2o = [2;0.7];
lambda2= eig(H_2i2o*H_2i2o');
C_2i2o = sum(log2(1+ lambda2*N0_inv/2));
C_2iso= log2(1+N0_inv/(2*(H_2iso*H_2iso')));
C_si2o = log2(1+N0_inv/(H_si2o'*H_si2o));
C_siso = log2(1+N0_inv/(H_siso*H_siso));
plot(SNR_db,[C_siso;C_2iso;C_si2o;C_2i2o]);
legend(['1x1';'2x1';'1x2';'2x2']);
xlabel('SNR');
ylabel('Capacity')

