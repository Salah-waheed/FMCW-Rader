% Task B: Range estametion (with CFAR for detection bouns)
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

% Setup Fig
fig_cfar = figure('Name', 'Task B: CFAR', 'Position', [50 50 1200 900], 'Color', 'w');

%% Main Loop
for i = 1:length(snr_vals)
    snr = snr_vals(i);
    
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
    
    % add noise
    srx_noisy = awgn(srx, snr, 'measured'); 
   
    %% Processing
    tic;
    
    % mixer and reshape
    mix = stx .* conj(srx_noisy);
    mix_mat = zeros(Npri, Np);
    for n = 1:Np
        k = (n-1)*Npri + 1;
        mix_mat(:,n) = mix(k : k+Npri-1).';
    end

    mix_mat = mix_mat .* hamming(Npri) .* hamming(Np)';
    
    %  window and range fft
    fft_r = fft(mix_mat, [], 1);
    fft_r = fft_r(1:Npri/2, :);
    fr = (0:Npri/2-1) * (Fs/Npri);
    r_axis = c * fr / (4 * slope); 
    
    %% CFAR
    cut = abs(fft_r(:,1));
    cut_db = 20*log10(cut);
    N_bins = length(cut_db);
    
    % settings
    Tr = 12;  
    Gd = 4;   
    os = 18; 
    
    Thresh = zeros(size(cut_db));
    
    for j = 1:N_bins
        p1 = max(1, j - (Tr + Gd));
        p2 = min(N_bins, j + (Tr + Gd));
        
        train_idx = [p1 : max(1, j - Gd - 1), min(N_bins, j + Gd + 1) : p2];
        
        if isempty(train_idx)
            noise = cut_db(j);
        else
            noise = mean(cut_db(train_idx));
        end
        
        Thresh(j) = noise + os;
    end
    
    % peak find
    [pks, locs] = findpeaks(cut_db, 'MinPeakProminence', 5);
    det_pks = [];
    det_locs = [];
    
    for k = 1:length(locs)
        idx = locs(k);
        if pks(k) > Thresh(idx)
            det_pks = [det_pks, pks(k)];
            det_locs = [det_locs, idx];
        end
    end
    
    det_r = r_axis(det_locs);
    t_proc = toc;
    
    % Output 
    fprintf('\n--- SNR %d dB ---\n', snr);
    disp('Targets Found:');
    for k = 1:nt
        r_true = tgt(k,1);
        if ~isempty(det_r)
            [err, ii] = min(abs(det_r - r_true));
            if err < 5
                fprintf(' Tgt %d: True=%.2f m, Det=%.2f m\n', k, r_true, det_r(ii));
            else
                fprintf(' Tgt %d: Not Detected\n', k);
            end
        else
            fprintf(' Tgt %d: Not Detected\n', k);
        end
    end
    fprintf(' Proc time: %.5f s\n', t_proc);
    
    %% Plots
    figure(fig_cfar);
    
    % Beat (Time)
    subplot(3, 2, (i-1)*2 + 1); 
    t_b = (0:Npri-1)*Ts*1e6;
    plot(t_b, real(mix_mat(:,1))); 
    title(sprintf('Beat (Time) SNR %d', snr));
    xlabel('Time (us)'); ylabel('Amp'); xlim([0 Tc*1e6]); grid on;
    
    % Range (Freq)
    subplot(3, 2, (i-1)*2 + 2);
    plot(r_axis, cut_db, 'b'); hold on;
    plot(r_axis, Thresh, 'g--', 'LineWidth', 1.5);
    
    if ~isempty(det_pks)
        plot(det_r, det_pks, 'ro', 'MarkerSize', 6, 'LineWidth', 1.5);
    end
    
    title(sprintf('Range (CFAR) SNR %d', snr));
    xlabel('Range (m)'); ylabel('dB'); xlim([0 250]); grid on;
    if i==1; legend('Sig', 'Thresh', 'Det'); end
end