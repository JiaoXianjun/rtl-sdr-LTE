% Jiao Xianjun (putaoshu@msn.com; putaoshu@gmail.com)
% Find out PSS location by moving FFT, peak averaging, Peak-to-Average-Ratio monitoring.
% A script of project: https://github.com/JiaoXianjun/rtl-sdr-LTE

function [hit_flag, hit_idx, hit_avg_snr, hit_snr, pss_idx] = move_fft_snr_runtime_avg(s, pss, th)
num_pss = size(pss, 2);

hit_flag = false;
hit_idx = -1;
pss_idx = -1;
hit_avg_snr = inf;
hit_snr = inf;

mv_len = 128+9; % the last OFDM symbol (SSS) of slot 1 and 11
distance_pss_sss = 9+10+9 + 3*128; % normal CP case
store_for_moving_avg = 999.*ones(num_pss, mv_len);
sampling_rate = 1.92e6; % LTE spec
fo_search_set = [-200e3:50e3:200e3];

pss_set = kron(ones(1, length(fo_search_set)), pss).*exp(1i.*2.*pi.*(1./sampling_rate).*(0:(length(pss)-1)).'*fo_search_set);

sum_snr = sum(store_for_moving_avg, 2);

len = length(s);
corr_len = size(pss, 1);
snr_record = zeros(num_pss, len - (corr_len-1));


signal_power = zeros(num_pss, 1);
for i=1:(len - (corr_len-1))
    chn_tmp = s(i:(i+corr_len-1));
    
    chn_tmp = pss'.*kron(ones(num_pss, 1), chn_tmp.');
    chn_tmp = abs(fft(chn_tmp, corr_len, 2)).^2;

%     tmp = fft(chn_tmp, corr_len);
%     tmp1 = [tmp; tmp(1:(end-1))];
%     tmp1_mat = lin2col_shift_mat(tmp1, corr_len);
%     tmp1_mat = [tmp1_mat(:, 1:(corr_len/4)) tmp1_mat(:, ((corr_len*3/4) + 1):end)]; % discard large frequency offset to save computations
%     chn_tmp = ( abs(fd_pss'*tmp1_mat).^2 );
    
%     signal_power = max(chn_tmp);
    [~, max_idx] = max(chn_tmp, [], 2);
    for j=1:num_pss
%         max_set = mod((max_idx(j) + (-1:1))-1, corr_len/2) + 1;
        max_set = mod((max_idx(j) + (-1:1))-1, corr_len) + 1;
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
%         chn_tmp = [chn_tmp(pss_idx, 1:(corr_len/4)), zeros(1, corr_len/2), chn_tmp(pss_idx, ((corr_len/4)+1) : end)];
        chn_tmp = chn_tmp(pss_idx,:);
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
