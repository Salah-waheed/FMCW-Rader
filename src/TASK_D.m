% Task D: range-doppler Map
% BY: Salah Waheed, Islam Ahmed
clear; close all; clc;

%% paramaters

c = 3e8; 
fc = 76.5e9; %carrier
lambda = c/fc;
B = 1e9; %Bandwidht
Tc = 3e-6; % Chirp duration
PRI = 8.4e-6; 
slope = B/Tc;
Ts = 0.5e-9; 
Fs = 1/Ts;
Np = 512;           % num pulses
Ns = round(Tc/Ts);   % samples per chirp
Npri = round(PRI/Ts); % samples per PRI

snr_vals = [5, 10, 15]; 

%% Targets [R, V]
tgt = [
    40    30;
    120  -80
];
nt = size(tgt,1);

% Axes
fr = (0:Npri/2-1) * (Fs/Npri);
r_ax = c * fr / (4 * slope); 
v_ax = (-Np/2 : Np/2-1) * (lambda / (2 * Np * PRI));

% Setup Fig
figure('Name', 'Task D: RD Maps', 'Position', [100 100 1200 400], 'Color', 'w');

%% Main Loop
for i = 1:length(snr_vals)
    snr = snr_vals(i);
    fprintf('Processing SNR %d dB...\n', snr);
    
    % Gen Signals
    N_all = Np * Npri;
    stx = zeros(1, N_all);
    srx = zeros(1, N_all);
    t_fast = (0:Ns-1)*Ts;
    
    % tx gen
    for n = 0:Np-1
        k = n*Npri + 1;
        stx(k:k+Ns-1) = exp(1j * pi * slope * t_fast.^2);
    end
    
    % rx gen
    for k = 1:nt
        R = tgt(k,1); 
        V = tgt(k,2); 
        Amp = (40 / R)^2; 
        
        for n = 0:Np-1
            t_slow = n * PRI;
            r_curr = R + V * t_slow;
            tau = 2 * r_curr / c;
            d_samp = round(tau/Ts);
            phs = exp(-1j * 2 * pi * fc * tau);
            
            idx = n*Npri + 1 + d_samp;
            
            if idx+Ns-1 <= N_all
                sig = Amp * (exp(1j * pi * slope * (t_fast - tau).^2) * phs);
                srx(idx:idx+Ns-1) = srx(idx:idx+Ns-1) + sig;
            end
        end
    end
    
    % noise
    srx_noisy = awgn(srx, snr, 'measured'); 
    
    %% Processing
    mix = stx .* conj(srx_noisy);
    
    % reshape and window
    mix_mat = zeros(Npri, Np);
    for n = 1:Np
        k = (n-1)*Npri + 1;
        mix_mat(:,n) = mix(k : k+Npri-1).';
    end

    mix_mat = mix_mat .* hamming(Npri) .* hamming(Np)';
    
    % 2D FFT
    fft_r = fft(mix_mat, [], 1);
    fft_r = fft_r(1:Npri/2, :); % range FFT
    
    fft_d = fftshift(fft(fft_r, [], 2), 2); % doppler FFT
    
    map_db = 20*log10(abs(fft_d) + eps); % Log scale
    map_norm = map_db - max(map_db(:));
    
    %% Plotting
    subplot(1, 3, i);
    imagesc(r_ax, v_ax, map_norm.'); 
    axis xy; colormap jet; 
    
    xlabel('Range (m)'); ylabel('Vel (m/s)');
    title(sprintf('SNR = %d dB', snr));
    clim([-75 0]); 
    
    hold on;
    plot(tgt(:,1), tgt(:,2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    
    if i == 1
        legend('Targets', 'Location', 'best', 'TextColor', 'w');
    end
end
