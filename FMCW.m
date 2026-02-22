% FMCW rader full mergad code with all tasks
% BY: Salah Waheed, Islam Ahmed
clear; close all; clc;

%% paramaters
c = 3e8; 
fc = 76.5e9; %carrier
lambda = c/fc;
B = 1e9; %Bandwidth
Tc = 3e-6; % Chirp duration
PRI = 8.4e-6; 
slope = B/Tc;
Ts = 0.5e-9; 
Fs = 1/Ts;
Np = 512;           % num pulses
Ns = round(Tc/Ts);   % samples per chirp
Npri = round(PRI/Ts); % samples per PRI

%% Targets [R, V]
tgt = [156, 10; 68, -80];
nt = size(tgt, 1);

%% instentous freq
figure;
n_show = 3;
t_ax = 0:Ts:(n_show*PRI - Ts);
tx_f = (fc - B/2) * ones(size(t_ax));

% tx line
for i = 0:n_show-1
    t0 = i*PRI;
    idx = (t_ax >= t0) & (t_ax < t0 + Tc);
    tx_f(idx) = (fc - B/2) + slope*(t_ax(idx) - t0);
end
plot(t_ax*1e6, tx_f/1e9, 'k', 'LineWidth', 2); hold on; grid on;

% rx lines
cols = lines(nt);
for i = 1:nt
    delay = 2 * tgt(i,1) / c;
    rx_f = (fc - B/2) * ones(size(t_ax));
    for j = 0:n_show-1
        t0 = j*PRI;
        idx = (t_ax >= t0 + delay) & (t_ax < t0 + delay + Tc);
        if any(idx)
            rx_f(idx) = (fc - B/2) + slope*(t_ax(idx) - t0 - delay);
        end
    end
    plot(t_ax*1e6, rx_f/1e9, 'Color', cols(i,:), 'LineWidth', 1.5);
end
xlabel('Time (us)'); ylabel('Freq (GHz)'); title('Freq vs Time');
xlim([0 n_show*PRI*1e6]);

%% Main Loop for snr values
snr_vals = [5, 10, 15];
for s = 1:length(snr_vals)
    snr = snr_vals(s);
    disp(['Processing SNR: ' num2str(snr)]);
    
    % sig setup
    N_all = Np * Npri;
    stx = zeros(1, N_all);
    srx = zeros(1, N_all);
    t_fast = (0:Ns-1)*Ts;
    
    % tx gen
    for i = 0:Np-1
        k = i*Npri + 1;
        stx(k:k+Ns-1) = exp(1j * pi * slope * t_fast.^2);
    end
    
    % rx gen
    for i = 1:nt
        R = tgt(i,1);
        V = tgt(i,2);
        Amp = (40/R)^2;
        
        for j = 0:Np-1
            t_slow = j * PRI;
            r_t = R + V*t_slow;
            tau = 2 * r_t / c;
            delay_n = round(tau/Ts);
            
            phs = exp(-1j*2*pi*fc*tau);
            k = j*Npri + 1 + delay_n;
            
            if k+Ns-1 <= N_all
                sig = Amp * (exp(1j*pi*slope*(t_fast - tau).^2) * phs);
                srx(k:k+Ns-1) = srx(k:k+Ns-1) + sig;
            end
        end
    end
    
    % add noise
    srx_noisy = awgn(srx, snr, 'measured');
    
    % Plot(TX/RX) for this SNR

    figure('Color','w', 'Name', ['Signals at SNR ' num2str(snr)]);
    time_plot = (0:Npri-1)*Ts*1e6; % microseconds
    
    subplot(2,1,1); 
    plot(time_plot, real(stx(1:Npri)), 'b'); 
    title('TX Signal 1 Chirp'); 
    xlabel('Time (us)'); ylabel('Amp'); axis tight; grid on;
    
    subplot(2,1,2); 
    plot(time_plot, real(srx_noisy(1:Npri)), 'r'); 
    title(['RX Signal at SNR ' num2str(snr) ' dB']); 
    xlabel('Time (us)'); ylabel('Amp'); axis tight; grid on;

    % Processing
    tic;
    mix = stx .* conj(srx_noisy);
    
    % reshape
    mix_mat = zeros(Npri, Np);
    for i = 1:Np
        k = (i-1)*Npri + 1;
        mix_mat(:,i) = mix(k : k+Npri-1).';
    end
    
    % window & fft
    win_2d = mix_mat .* hamming(Npri) .* hamming(Np)';
    
    fft_r = fft(win_2d, [], 1);
    fft_r = fft_r(1:Npri/2, :); % one sided
    
    % axes
    fr = (0:Npri/2-1)*(Fs/Npri);
    r_axis = c * fr / (4 * slope);
    
    %% CFAR
    % detection on 1st chirp
    cut = abs(fft_r(:,1));
    cut_db = 20*log10(cut);
    N_bins = length(cut_db);
    
    % settings
    Tr = 12; Gd = 5; off = 17;
    Thresh = zeros(size(cut_db));
    
    for i = 1:N_bins
        i1 = max(1, i-(Tr+Gd));
        i2 = min(N_bins, i+(Tr+Gd));
        
        % noise window
        win_idx = [i1 : max(1, i-Gd-1), min(N_bins, i+Gd+1) : i2];
        
        noise = mean(cut_db(win_idx));
        Thresh(i) = noise + off;
    end
    
    % find peaks
    [pks, locs] = findpeaks(cut_db, 'MinPeakProminence', 5);
    det_pks = []; 
    det_locs = [];
    
    for i = 1:length(locs)
        idx = locs(i);
        if pks(i) > Thresh(idx)
            det_pks = [det_pks, pks(i)];
            det_locs = [det_locs, idx];
        end
    end
    det_r = r_axis(det_locs);
    t_end = toc;
    
    %% Plotting Results
    figure('Color','w', 'Name', ['Results SNR ' num2str(snr)]);
    
    % beat sig
    subplot(2,2,1);
    t_b = (0:Npri-1)*Ts;
    plot(t_b*1e6, real(mix_mat(:,1))); axis tight;
    xlabel('Time (us)'); ylabel('Amp'); title('Beat Signal');
    xlim([0 Tc*1e6]);
    
    % range-mag
    subplot(2,2,2); hold on; grid on;
    plot(r_axis, cut_db, 'b');
    plot(r_axis, Thresh, 'g--');
    if ~isempty(det_r)
        plot(det_r, det_pks, 'ro', 'LineWidth', 2);
    end
    xlabel('Range (m)'); ylabel('dB'); title(['Range at SNR ' num2str(snr)]);
    xlim([0 250]);
    
    % print results
    disp('Targets Found:');
    for i = 1:nt
        r_true = tgt(i,1);
        if ~isempty(det_r)
            [err, idx] = min(abs(det_r - r_true));
            if err < 5
                fprintf(' Tgt %d: True=%.2f m, Det=%.2f m\n', i, r_true, det_r(idx));
            else
                fprintf(' Tgt %d: Not Detected\n', i);
            end
        else
            fprintf(' Tgt %d: Not Detected\n', i);
        end
    end
    fprintf(' Proc time: %.5f s\n', t_end);
    
    % velocity-mag
    v_axis = (-Np/2 : Np/2-1) * (lambda / (2*Np*PRI));
    n_det = length(det_r);
    
    if n_det > 0
        subplot(2,2,3); hold on; grid on;
        disp('Velocity:');
        for i = 1:n_det
            idx = det_locs(i);
            % fft on slow time
            s_dopp = fft_r(idx, :);
            f_dopp = fftshift(fft(s_dopp));
            p_dopp = 20*log10(abs(f_dopp));
            
            [~, imax] = max(p_dopp);
            v_meas = v_axis(imax);
            
            plot(v_axis, p_dopp);
            plot(v_meas, p_dopp(imax), 'ro');
            
            % print velocety
            curr_r = det_r(i);
            [~, tidx] = min(abs(tgt(:,1) - curr_r));
            v_true = tgt(tidx, 2);
            fprintf(' R=%.1f m: True=%.2f m/s, Meas=%.2f m/s\n', curr_r, v_true, v_meas);
        end
        xlabel('m/s'); ylabel('dB'); title('velocity');
    else
        subplot(2,2,3); title('No Velocity');
    end
    
    % 2D map
    fft2 = fftshift(fft(fft_r, [], 2), 2);
    map = 20*log10(abs(fft2) + eps);
    map = map - max(map(:));
    
    subplot(2,2,4); 
    imagesc(r_axis, v_axis, map.'); 
    axis xy; colormap jet; clim([-70 0]);
    xlabel('Range'); ylabel('Vel'); title('Range-velocity Map');
    hold on;
    plot(tgt(:,1), tgt(:,2), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    colorbar;
end