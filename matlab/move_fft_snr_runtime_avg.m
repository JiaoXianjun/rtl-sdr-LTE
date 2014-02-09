% Jiao Xianjun (putaoshu@msn.com; putaoshu@gmail.com)
% Find out PSS location by moving FFT, peak averaging, Peak-to-Average-Ratio monitoring.
% A script of project: https://github.com/JiaoXianjun/rtl-sdr-LTE

function [hit_flag, hit_idx, hit_avg_snr, hit_snr, hit_fo, pss_idx] = move_fft_snr_runtime_avg(s, fd_pss, mv_len, fft_len, th)
num_pss = size(fd_pss, 2);

hit_flag = false;
hit_idx = -1;
pss_idx = -1;
hit_avg_snr = inf;
hit_snr = inf;
hit_fo = inf;

store_for_moving_avg = 999.*ones(num_pss, mv_len);
sum_snr = sum(store_for_moving_avg, 2);

len = length(s);
snr_record = zeros(num_pss, len - (fft_len-1));

sampling_rate = 1.92e6; % LTE spec

signal_power = zeros(num_pss, 1);
for i=1:(len - (fft_len-1))
    chn_tmp = s(i:(i+fft_len-1));
    
%     chn_tmp = abs(fft(chn_tmp, fft_len)).^2;
    tmp = fft(chn_tmp, fft_len);
    tmp1 = [tmp; tmp(1:(end-1))];
    tmp1_mat = lin2col_shift_mat(tmp1, fft_len);
    tmp1_mat = [tmp1_mat(:, 1:(fft_len/4)) tmp1_mat(:, ((fft_len*3/4) + 1):end)]; % discard large frequency offset to save computations
    chn_tmp = ( abs(fd_pss'*tmp1_mat).^2 );
    
%     signal_power = max(chn_tmp);
    [~, max_idx] = max(chn_tmp, [], 2);
    for j=1:num_pss
        max_set = mod((max_idx(j) + (-1:1))-1, fft_len/2) + 1;
%         max_set = mod((max_idx(j) + (-1:1))-1, fft_len) + 1;
        signal_power(j) = sum( chn_tmp(j, max_set) );
    end

    noise_power = sum(chn_tmp, 2) - signal_power;
    snr = 10.*log10(signal_power./noise_power);
    snr_record(:, i) = snr;
    
    peak_to_avg = snr - (sum_snr./mv_len);

    logic_tmp = (peak_to_avg > th);
    if sum(logic_tmp)
        hit_flag = true;
        pss_idx = find(logic_tmp==1);
        chn_tmp = [chn_tmp(pss_idx, 1:(fft_len/4)), zeros(1, fft_len/2), chn_tmp(pss_idx, ((fft_len/4)+1) : end)];
%         chn_tmp = chn_tmp(pss_idx,:);
        hit_fo = fo_from_fft_result(chn_tmp.', sampling_rate);
        hit_idx = i;
        hit_snr = snr(pss_idx);
        hit_avg_snr = snr(pss_idx) - peak_to_avg(pss_idx);
        figure; plot(chn_tmp);
%         disp(['Hit. idx ' num2str(i) '; SNR ' num2str(snr) 'dB; peak SNR to avg SNR ' num2str(peak_to_avg) 'dB']);
        break;
    else
        sum_snr = sum_snr - store_for_moving_avg(:, end);
        sum_snr = sum_snr + snr;
        
        store_for_moving_avg(:, 2:end) = store_for_moving_avg(:, 1:(end-1));
        store_for_moving_avg(:, 1) = snr;
    end
    
end

figure; plot(snr_record.');
