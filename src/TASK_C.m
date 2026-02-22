% Task C: Velocity estametion ( with findpeaks() for detection bouns)
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

% Velocity Axis
v_ax = (-Np/2 : Np/2-1) * (lambda / (2 * Np * PRI));

%% Main Loop

for i = 1:length(snr_vals)
    snr = snr_vals(i);
    
    fprintf('\n--- SNR %d dB ---\n', snr);
    
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
    
    % range fft
    fft_r = fft(mix_mat, [], 1);
    fft_r = fft_r(1:Npri/2, :);
    fr = (0:Npri/2-1) * (Fs/Npri);
    r_ax = c * fr / (4 * slope); 
    
    % Peak Detection 
    mag_db = 20*log10(abs(fft_r(:,1)));
    thresh = max(mag_db) - 25; 
    [~, locs] = findpeaks(mag_db, 'MinPeakHeight', thresh, 'MinPeakProminence', 5);
    det_r = r_ax(locs);
    n_det = length(det_r);
    figure('Name', sprintf('SNR %d Velocity', snr), 'Color', 'w');
    
    % Loop detected targets & Print Results
    if n_det > 0
        disp('Targets Found:');
        for j = 1:n_det
            % Extract bin & FFT
            r_idx = locs(j);
            slow_sig = fft_r(r_idx, :); 
            spec = fftshift(fft(slow_sig));
            spec_db = 20*log10(abs(spec));
            
            % Find Peak Velocity
            [~, pks_loc] = findpeaks(spec_db, 'SortStr', 'descend', 'NPeaks', 1);
            if isempty(pks_loc)
                [~, max_i] = max(spec_db); 
            else
                max_i = pks_loc(1); 
            end
            v_est = v_ax(max_i);
            
            % Match to truth
            r_curr = det_r(j);
            [~, ti] = min(abs(tgt(:,1) - r_curr));
            v_true = tgt(ti, 2);
            fprintf(' Range=%.1fm: True=%.2f m/s, Det=%.2f m/s\n', r_curr, v_true, v_est);
            subplot(n_det, 1, j);
            plot(v_ax, spec_db, 'b', 'LineWidth', 1.5); hold on;
            plot(v_est, spec_db(max_i), 'ro', 'MarkerSize', 6, 'LineWidth', 2);
            xline(v_est, 'r--');          
            title(sprintf('Target @ %.1fm', r_curr));
            xlabel('Vel (m/s)'); ylabel('dB'); grid on; xlim([-100 100]);
        end
    else
        disp('No Targets Detected');
    end
end