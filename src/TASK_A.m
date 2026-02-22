% Task A: TX,RX GEN
% BY: Salah Waheed, Islam Ahmed
clear; close all; clc;

%% Paramaters

c = 3e8; 
fc = 76.5e9;   %carrier       
B = 1e9; %Bandwidht
Tc = 3e-6; % chirp duration
PRI = 8.4e-6;
slope = B/Tc;
f_start = fc - B/2; 
f_stop = fc + B/2;
Ts = 0.5e-9;
Fs = 1/Ts;
Np = 512;           % num pulses
Ns = round(Tc/Ts);   % samples per chirp
Npri = round(PRI/Ts); % samples per PRI
snr_vals = [5, 10, 15]; 

%% Targets [R, V]
Tgt = [
    40    30;
    120  -80
];
n_tgt = size(Tgt,1);

%% instentous freq
figure('Name', 'Freq Check', 'Position', [50 50 800 400], 'Color', 'w');
n_show = 3;
t_ax = 0:Ts:(n_show*PRI - Ts);
len_t = length(t_ax);

% Tx line
tx_freq = f_start * ones(size(t_ax));
for n = 0:n_show-1
    t0 = n*PRI;
    mask = (t_ax >= t0) & (t_ax < t0 + Tc);
    tx_freq(mask) = f_start + slope*(t_ax(mask) - t0);
end

% Rx lines
rx_freq = f_start * ones(n_tgt, len_t); 
for k = 1:n_tgt
    delay = 2 * Tgt(k,1) / c;
    for n = 0:n_show-1
        t0 = n*PRI;
        mask = (t_ax >= t0 + delay) & (t_ax < t0 + delay + Tc);
        if any(mask)
            rx_freq(k, mask) = f_start + slope*(t_ax(mask) - t0 - delay);
        end
    end
end

% Plotting
plot(t_ax*1e6, tx_freq/1e9, 'k', 'LineWidth', 2); hold on;
cols = lines(n_tgt);
lgd = {'TX'};

for k = 1:n_tgt
    plot(t_ax*1e6, rx_freq(k,:)/1e9, 'LineWidth', 1.5, 'Color', cols(k,:));
    lgd{end+1} = sprintf('RX Tgt %d', k);
end

xlabel('Time (us)'); ylabel('Freq (GHz)');
title('Instantaneous Freq');
legend(lgd, 'Location', 'bestoutside');
ylim([75.8 77.2]); xlim([0 n_show*PRI*1e6]); grid on;

%% TX,RX signals 
% Figure setup
fig_raw = figure('Name', 'Raw Signals', 'Color', 'w', 'Position', [100 100 800 900]);
t_sec = (0:Npri-1)*Ts;

% Generate TX 

Nt = Np * Npri;
Tx = zeros(1, Nt);
t_base = (0:Ns-1)*Ts;

for n = 0:Np-1
    idx = n*Npri + 1;
    Tx(idx:idx+Ns-1) = exp(1j * pi * slope * t_base.^2);
end

figure(fig_raw);
subplot(4,1,1); 
plot(t_sec*1e6, real(Tx(1:Npri)), 'b');
title('TX');
xlabel('Time (us)'); ylabel('Amp');
grid on; axis tight;

% Loop RX for SNRs
for i = 1:length(snr_vals)
    snr = snr_vals(i);
    fprintf('Running SNR %d dB...\n', snr);

    % Rx gen
    Rx = zeros(1, Nt);
    for k = 1:n_tgt
        R = Tgt(k,1); 
        V = Tgt(k,2); 
        A = (40 / R)^2; 
        
        for n = 0:Np-1
            t_slow = n * PRI;
            r_curr = R + V * t_slow;
            tau = 2 * r_curr / c;
            d_samp = round(tau/Ts);
            phs = exp(-1j * 2 * pi * fc * tau);
            
            idx = n*Npri + 1 + d_samp;
            
            if idx+Ns-1 <= Nt
                sig = A * (exp(1j * pi * slope * (t_base - tau).^2) * phs);
                Rx(idx:idx+Ns-1) = Rx(idx:idx+Ns-1) + sig;
            end
        end
    end

    % Add noise
    Rx_noisy = awgn(Rx, snr, 'measured'); 

    % Plotting
    figure(fig_raw);
    subplot(4,1,i+1); 
    plot(t_sec*1e6, real(Rx_noisy(1:Npri)), 'r');
    title(sprintf('RX at SNR %d dB', snr));
    xlabel('Time (us)'); ylabel('Amp');
    grid on; axis tight;
end