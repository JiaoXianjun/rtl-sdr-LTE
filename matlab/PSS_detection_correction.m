% Jiao Xianjun (putaoshu@msn.com; putaoshu@gmail.com)
% Find out LTE PSS in the signal stream and correct sampling&carrier error.
% A script of project: https://github.com/JiaoXianjun/rtl-sdr-LTE

function [position, pss_idx] = PSS_detection_correction(s, fd_pss, td_pss)
disp(' ');
position = -1;
pss_idx = -1;

fft_len = size(td_pss, 1);

len = length(s);
th = 100; %dB. threshold

sampling_rate = 1.92e6; % LTE spec. 30.72MHz/16.

len_time_subframe = 1e-3; % 1ms. LTE spec
num_sample_per_subframe = len_time_subframe*sampling_rate;
num_subframe_per_radioframe = 10;
num_sample_per_radioframe = num_sample_per_subframe*num_subframe_per_radioframe;

% find out first PSS in the first radio frame by moving FFT
[hit_flag, hit_idx, hit_avg_snr, hit_snr, pss_idx] = move_fft_snr_runtime_avg(s(1:num_sample_per_radioframe), td_pss((end-fft_len+1):end,:), th);
disp([hit_idx hit_fo]);
if ~hit_flag
    disp('PSS coarse: No PSS found!');
    return;
end

num_sym_between_FCCH = 10*num_slot_per_frame*num_sym_per_slot;
num_sym_between_FCCH1 = 11*num_slot_per_frame*num_sym_per_slot; % in case the last idle frame of the multiframe

num_sym_between_FCCH_decimate = round(num_sym_between_FCCH/decimation_ratio);
num_sym_between_FCCH_decimate1 = round(num_sym_between_FCCH1/decimation_ratio);

max_num_fcch = ceil(len/(10*num_sym_per_frame/decimation_ratio));
position = zeros(1, max_num_fcch);
snr = zeros(1, max_num_fcch);
position(1) = hit_idx;
snr(1) = hit_snr;

set_idx = 1;
max_offset = 5;
while 1
    next_position = position(set_idx) + num_sym_between_FCCH_decimate; % predicted position of next FCCH in the same multiframe
    
    if next_position > (len - (fft_len-1)) - max_offset; % run out of sampled signal
        break;
    end

%     i_set = next_position + (-max_offset:max_offset);
%     i_set(i_set<1) = 1;
%     i_set(i_set>(len - (fft_len-1))) = (len - (fft_len-1));
    i_set = next_position + [-max_offset, max_offset];

    [hit_flag, hit_idx, hit_snr] = specific_fft_snr_fix_avg(s, i_set, fft_len, th, hit_avg_snr); % fft detection at specific position
    
    if hit_flag
        set_idx = set_idx + 1;
        position(set_idx) = hit_idx;
        snr(set_idx) = hit_snr;
    else
        next_position = position(set_idx) + num_sym_between_FCCH_decimate1;% predicted position of next FCCH in the next multiframe
        
        if next_position > (len - (fft_len-1)) - max_offset; % run out of sampled signal
            break;
        end

%         i_set = next_position + (-max_offset:max_offset);
%         i_set(i_set<1) = 1;
%         i_set(i_set>(len - (fft_len-1))) = (len - (fft_len-1));
        i_set = next_position + [-max_offset, max_offset];
        
        [hit_flag, hit_idx, hit_snr] = specific_fft_snr_fix_avg(s, i_set, fft_len, th, hit_avg_snr); % fft detection at specific position
        
        if hit_flag
            set_idx = set_idx + 1;
            position(set_idx) = hit_idx;
            snr(set_idx) = hit_snr;
        else
            break;
        end
    end
end

position = position(1:set_idx);
snr = snr(1:set_idx);

position = (position-1)*decimation_ratio + 1;
disp(['PSS coarse: hit successive ' num2str(length(position)) ' FCCH. pos ' num2str(position)]);
disp(['PSS coarse: pos diff ' num2str(diff(position))]);
disp(['PSS coarse: SNR ' num2str(snr)]);
