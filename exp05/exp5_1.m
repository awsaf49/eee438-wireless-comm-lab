clc;
clear all;
close all;

consts = [1+1i, 1-1i, -1+1i, -1-1i];
Es = 1;

SNR_db = -10:5:20;
N0 = Es./10.^(SNR_db/10);
sim_len = 1e5;
num_sym = sim_len*2;

figure(1);

SER = zeros(1, length(SNR_db));
for idx=1:length(SNR_db)
    num_errs = 0;
    for sim_no = 1:sim_len
        s_tx = 2*randi([0, 1], 2, 1)-1 + 1i*(2*randi([0, 1], 2, 1)-1);
        H = 1/sqrt(2) * (randn(1, 2) + 1i * randn(1, 2));
        noise = sqrt(N0(idx))* (randn(2, 1) + 1i*randn(2, 1));
        H22 = zeros(2,2);
        H22(1, 1) = H(1, 1);
        H22(1, 2) = H(1, 2);
        H22(2, 1) = H(1, 2);
        H22(2, 2) = -H(1, 1);
        y =  H22 * s_tx + noise;
        s_est = H22' * y;
        dist1 = abs(s_est(1) - consts);
        dist2 = abs(s_est(2) - consts);
        [~, idx1] = min(dist1);
        [~, idx2] = min(dist2);
        s_rx(1, 1) = consts(idx1);
        s_rx(2, 1) = consts(idx2);
        errs = (s_rx ~= s_tx);
        num_errs = num_errs + sum(errs(:));
        
    end
    SER(idx) = num_errs / num_sym;
end
semilogy(SNR_db, SER);
hold on;

%% part2
SER = zeros(1, length(SNR_db));
for idx=1:length(SNR_db)
    num_errs = 0;
    for sim_no = 1:sim_len
        s_tx = 2*randi([0, 1], 2, 1)-1 + 1i*(2*randi([0, 1], 2, 1)-1);
        H = 1/sqrt(2) * (randn(2, 2) + 1i * randn(2, 2));
        noise = sqrt(N0(idx))* (randn(4, 1) + 1i*randn(4, 1));
        H22 = zeros(4,2);
        H22(1, 1) = H(1, 1);
        H22(1, 2) = H(1, 2);
        H22(2, 1) = H(2, 1);
        H22(2, 2) = H(2, 2);
        H22(3, 1) = H(1, 2);
        H22(3, 2) = -H(1, 1);
        H22(4, 1) = H(2, 2);
        H22(4, 2) = -H(2, 1);
        y =  H22 * s_tx + noise;
        s_est = H22' * y;
        dist1 = abs(s_est(1) - consts);
        dist2 = abs(s_est(2) - consts);
        [~, idx1] = min(dist1);
        [~, idx2] = min(dist2);
        s_rx(1, 1) = consts(idx1);
        s_rx(2, 1) = consts(idx2);
        errs = (s_rx ~= s_tx);
        num_errs = num_errs + sum(errs(:));
        
    end
    SER(idx) = num_errs / num_sym;
end
semilogy(SNR_db, SER);
hold on;
grid on;
xlabel('SNR(dB)');
ylabel('SER');
title('SER vs SNR');
legend('2x1', '2x2');
