% Jiao Xianjun (putaoshu@msn.com; putaoshu@gmail.com)
% Find out LTE PSS in the signal stream and correct sampling&carrier error.
% A script of project: https://github.com/JiaoXianjun/rtl-sdr-LTE

function [ppm, f_set] = sampling_ppm_f_search_set_by_pss(s, td_pss)
% sampling period PPM! not sampling frequency PPM!

len_pss = size(td_pss, 1);

ppm = inf;
f_set = inf;

fo_search_set = -100e3 : 5e3 : 100e3; % -100kHz ~ 100 kHz with 5kHz step size
pss_fo_set = pss_fo_set_gen(td_pss, fo_search_set);

len = length(s);
th = 25; %threshold

sampling_rate = 1.92e6; % LTE spec. 30.72MHz/16.

len_time_subframe = 1e-3; % 1ms. LTE spec
num_subframe_per_radioframe = 10;
num_sample_per_subframe = len_time_subframe*sampling_rate;
num_sample_per_radioframe = num_sample_per_subframe*num_subframe_per_radioframe;

[hit_pss_fo_set_idx, hit_time_idx, corr_val] = pss_moving_corr(s(1:4*num_sample_per_radioframe), pss_fo_set, th);
% disp(num2str(hit_pss_fo_set_idx));
% disp(num2str(hit_time_idx));

if ~sum(hit_pss_fo_set_idx)
    disp('No strong enough PSS correlation peak.');
    return;
end

pss_period = num_sample_per_radioframe/2;

num_fo = length(hit_pss_fo_set_idx);
max_num_hit = ceil(len/pss_period);
time_location = zeros(max_num_hit, num_fo);
time_location(1,:) = hit_time_idx;
% time_location_invalid_record = zeros(max_num_hit, num_fo);
hit_corr_val = zeros(max_num_hit, num_fo);
hit_corr_val(1,:) = corr_val;

pss_count = 1;
max_offset = 32;
while 1
    next_location = time_location(pss_count,:) + pss_period; % predicted position of next PSS
    min_next_location = min(next_location);
    max_next_location = max(next_location);
%     disp(num2str([max_next_location-min_next_location min_next_location max_next_location]));
    
    if max_next_location + max_offset > (len - (len_pss-1)); % run out of sampled signal
        break;
    end

    i_set = [min_next_location-max_offset, max_next_location+max_offset];

    [hit_time_idx, corr_val] = pss_fix_location_corr(s, i_set, pss_fo_set, hit_pss_fo_set_idx);
    
    pss_count = pss_count + 1;
    time_location(pss_count,:) = hit_time_idx;
    hit_corr_val(pss_count,:) = corr_val;
    
%     invalid_idx = find(corr_val<(max(corr_val)/2));
%     invalid_idx = find( (corr_val - (hit_corr_val(pss_count-1,:).*1./2) )<0 );
%     time_location(pss_count, invalid_idx) = next_location(invalid_idx);
    
%     time_location_invalid_record(pss_count, invalid_idx) = 1;
%     disp(['invalid idx ' num2str(invalid_idx)]);
end

time_location = time_location(1:pss_count, :);
hit_corr_val = hit_corr_val(1:pss_count,:);

diff_time_location = diff(time_location, 1, 1);
diff_time_location = abs( diff_time_location - pss_period );
[first_valid_idx, last_valid_idx] = find_best_timing_location(diff_time_location, round( (len/pss_period)*(1/2) ));
if first_valid_idx == -1
    disp('No long enough continuous PSS hit.');
    return;
end

ppm_store = zeros(1, num_fo);
fo_idx = find(first_valid_idx>0);
valid_idx = zeros(1, num_fo);
ppm_idx = 0;
for i=1:length(fo_idx)
% for i=1:1
    col_idx = fo_idx(i);
    sp = first_valid_idx(col_idx);
    ep = last_valid_idx(col_idx);
    
%     tmp_invalid_record = time_location_invalid_record(sp:ep, col_idx);
%     if sum(tmp_invalid_record) > 0
%         continue;
%     end
    
    tmp_time_idx = time_location(sp:ep, col_idx);
    
    distance = tmp_time_idx(end) - tmp_time_idx(1);
    len = length(tmp_time_idx);
    ppm = 1e6*( distance - (pss_period*(len-1)) )./(pss_period*(len-1)); % sampling period PPM! not sampling frequency PPM!

    ppm_idx = ppm_idx + 1;
    ppm_store(ppm_idx) = ppm;
    valid_idx(ppm_idx) = i;
end
if ppm_idx == 0
    disp('No valid PSS hit sequence.');
    return;
end
ppm = mean(ppm_store(1:ppm_idx));
disp(['Total ' num2str(ppm_idx) ' freq. idx for PPM: ' num2str(fo_idx(valid_idx(1:ppm_idx)))]);

num_f_reserve = 1;

sum_corr_val = sum(hit_corr_val, 1);
% disp(['Pre-Proc: sum corr ' num2str(sum_corr_val)]);

num_reserve = min([num_f_reserve, num_fo]);
[~, max_idx] = sort(sum_corr_val, 'descend');
max_idx = max_idx(1:num_reserve);
disp(['Freq. idx for f_set: ' num2str(max_idx)]);

idx_in_fo_search_set = hit_pss_fo_set_idx(max_idx);
f_set = fo_search_set(mod(idx_in_fo_search_set-1, length(fo_search_set)) + 1);

disp(['Period PPM ' num2str(ppm) 'PPM; f_set ' num2str(f_set./1e3) 'kHz']);

% % ----------------------------------------------------------------------

% time_location = time_location(1:pss_count, :);
% hit_corr_val = hit_corr_val(1:pss_count, :);
% % disp(num2str(diff(time_location, 1, 1)));
% 
% sum_corr_val = sum(hit_corr_val, 1);
% [~, max_idx] = max(sum_corr_val);
% sp = time_location(1, max_idx);
% ep = time_location(end, max_idx);
% ppm = 1e6*( (ep-sp) - (pss_period*(pss_count-1)) )./(pss_period*(pss_count-1)); % sampling period PPM! not sampling frequency PPM!
% 
% num_reserve = min(4, num_fo);
% [~, max_idx] = sort(sum_corr_val, 'descend');
% max_idx = max_idx(1:num_reserve);
% idx_in_fo_search_set = hit_pss_fo_set_idx(max_idx);
% f_set = fo_search_set(mod(idx_in_fo_search_set-1, length(fo_search_set)) + 1);
% disp(['Pre-Proc: Period PPM ' num2str(ppm) 'PPM; f_set ' num2str(f_set./1e3) 'kHz']);

% % ----------------------------------------------------------------------

% num_f_reserve = 4;
% sum_corr_val = sum(hit_corr_val, 1);
% [~, max_idx] = max(sum_corr_val);
% sp = time_location(1, max_idx);
% ep = time_location(end, max_idx);
% ppm = 1e6*( (ep-sp) - (pss_period*(pss_count-1)) )./(pss_period*(pss_count-1)); % sampling period PPM! not sampling frequency PPM!
% 
% num_reserve = min(num_f_reserve, num_fo);
% [~, max_idx] = sort(sum_corr_val, 'descend');
% max_idx = max_idx(1:num_reserve);
% idx_in_fo_search_set = hit_pss_fo_set_idx(max_idx);
% f_set = fo_search_set(mod(idx_in_fo_search_set-1, length(fo_search_set)) + 1);
% disp(['Pre-Proc: Period PPM ' num2str(ppm) 'PPM; f_set ' num2str(f_set./1e3) 'kHz']);
