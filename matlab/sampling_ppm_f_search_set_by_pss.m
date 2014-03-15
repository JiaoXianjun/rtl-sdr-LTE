% Jiao Xianjun (putaoshu@msn.com; putaoshu@gmail.com)
% Find out LTE PSS in the signal stream and correct sampling&carrier error.
% A script of project: https://github.com/JiaoXianjun/rtl-sdr-LTE

function [ppm, f_set] = sampling_ppm_f_search_set_by_pss(s, td_pss)
% sampling period PPM! not sampling frequency PPM!

len_pss = size(td_pss, 1);
num_pss = 3;

ppm = inf;
f_set = inf;

fo_search_set = -100e3 : 5e3 : 100e3; % -100kHz ~ 100 kHz with 5kHz step size
pss_fo_set = pss_fo_set_gen(td_pss, fo_search_set);

len = length(s);
th = 25*265.1154; %threshold

sampling_rate = 1.92e6; % LTE spec. 30.72MHz/16.

len_time_subframe = 1e-3; % 1ms. LTE spec
num_subframe_per_radioframe = 10;
num_sample_per_subframe = len_time_subframe*sampling_rate;
num_sample_per_radioframe = num_sample_per_subframe*num_subframe_per_radioframe;

[hit_pss_fo_set_idx, hit_time_idx, corr_val] = pss_moving_corr(s(1:(2*num_sample_per_radioframe)), pss_fo_set, th);
% disp(num2str(hit_pss_fo_set_idx));
% disp(num2str(hit_time_idx));

if ~sum(hit_pss_fo_set_idx)
    disp('No strong enough PSS correlation peak.');
    return;
end

num_fo = length(hit_pss_fo_set_idx);
max_reserve_per_pss = 8;
pss_idx = floor( (hit_pss_fo_set_idx-1)./length(fo_search_set) ) + 1;

% method for C
pss_reserve_idx_bin = zeros(num_pss,num_fo);
for i=1:num_pss
    pss_reserve_idx_bin(i,:) = (pss_idx==i);
    tmp_num = sum(pss_reserve_idx_bin(i,:));
    if tmp_num>max_reserve_per_pss
        num_discard = tmp_num - max_reserve_per_pss;
        for j=num_fo:-1:1
            if pss_reserve_idx_bin(i,j)==1
                pss_reserve_idx_bin(i,j) = 0;
                num_discard = num_discard - 1;
                if num_discard == 0
                    break;
                end
            end
        end
    end
end
reserve_idx = find(sum(pss_reserve_idx_bin, 1));

hit_pss_fo_set_idx = hit_pss_fo_set_idx(reserve_idx);
hit_time_idx = hit_time_idx(reserve_idx);
corr_val = corr_val(reserve_idx);

num_fo = length(hit_pss_fo_set_idx);

pss_period = num_sample_per_radioframe/2;

max_num_hit = ceil(len/pss_period);
time_location = zeros(max_num_hit, num_fo);
time_location(1,:) = hit_time_idx;
hit_corr_val = zeros(max_num_hit, num_fo);
hit_corr_val(1,:) = corr_val;

pss_count = 1;
max_offset = 32;
time_location_invalid_record = zeros(max_num_hit, num_fo);
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
%     plot(dbgtmp);
    
    pss_count = pss_count + 1;
    time_location(pss_count,:) = hit_time_idx;
    hit_corr_val(pss_count,:) = corr_val;
    
    invalid_logic = (corr_val<(th*3/4));
    time_location(pss_count, invalid_logic) = next_location(invalid_logic);
    
    time_location_invalid_record(pss_count, :) = invalid_logic;
%     disp(['invalid idx ' num2str(invalid_idx)]);
end

time_location = time_location(1:pss_count, :);
hit_corr_val = hit_corr_val(1:pss_count,:);
time_location_invalid_record = time_location_invalid_record(1:pss_count,:);

% diff_time_location = diff(time_location, 1, 1);
% diff_time_location = abs( diff_time_location - pss_period );
% [first_valid_idx, last_valid_idx] = find_best_timing_location(diff_time_location, round( (len/pss_period)*(1/2) ));
% if first_valid_idx == -1
%     disp('No long enough continuous PSS hit.');
%     return;
% end

ppm_store = zeros(1, num_fo);
% fo_idx = find(first_valid_idx>0);
valid_idx = zeros(1, num_fo);
min_dist = floor( (len/pss_period)*(1/2) );
ppm_idx = 0;
for i=1:num_fo
% for i=1:1
    col_idx = i;
    sp = find(time_location_invalid_record(:,col_idx)==0, 1, 'first');
    ep = find(time_location_invalid_record(:,col_idx)==0, 1, 'last');
    
%     tmp_invalid_record = time_location_invalid_record(sp:ep, col_idx);
%     if sum(tmp_invalid_record) > 0
%         continue;
%     end
    
    if isempty(sp)
        continue;
    end
    
    if (ep-sp)<min_dist
%     if  sum(time_location_invalid_record(:,col_idx)==0)<min_dist
        continue;
    end
    
    distance = time_location(ep,col_idx) - time_location(sp,col_idx);
    len = ep-sp+1;
    ppm = 1e6*( distance - (pss_period*(len-1)) )./(pss_period*(len-1)); % sampling period PPM! not sampling frequency PPM!

    ppm_idx = ppm_idx + 1;
    ppm_store(ppm_idx) = ppm;
    valid_idx(ppm_idx) = i;
end
if ppm_idx == 0
    disp('No valid PSS hit sequence.');
    return;
end
ppm_store = ppm_store(1:ppm_idx);
valid_idx = valid_idx(1:ppm_idx);

valid_idx_backup = valid_idx;

extra_frequency_flag = 0;
disp(['PPM: ' num2str(ppm_store)]);
pss_idx = floor( (hit_pss_fo_set_idx(valid_idx)-1)./length(fo_search_set) ) + 1;
disp(['PSS: ' num2str(pss_idx)]);
if length(ppm_store) == 1
    ppm = ppm_store;
    disp(['Total ' num2str(ppm_idx) ' freq. idx for PPM: ' num2str(valid_idx)]);
    pss_idx = floor( (hit_pss_fo_set_idx(valid_idx)-1)./length(fo_search_set) ) + 1;
    disp(['Total ' num2str(ppm_idx) '  pss. idx for PPM: ' num2str(pss_idx)]);
    disp(['Average PPM: ' num2str(ppm)]);
    idx_in_fo_search_set = hit_pss_fo_set_idx(valid_idx);
    f_set = fo_search_set(mod(idx_in_fo_search_set-1, length(fo_search_set)) + 1);
    disp(['Period PPM ' num2str(ppm) 'PPM; f_set ' num2str(f_set./1e3) 'kHz']);
    disp(['Final PSS idx ' num2str(pss_idx)]);
    return;
elseif length(ppm_store) == 2
    ppm = mean(ppm_store);
    disp(['Total ' num2str(ppm_idx) ' freq. idx for PPM: ' num2str(valid_idx)]);
    pss_idx = floor( (hit_pss_fo_set_idx(valid_idx)-1)./length(fo_search_set) ) + 1;
    disp(['Total ' num2str(ppm_idx) '  pss. idx for PPM: ' num2str(pss_idx)]);
    disp(['Average PPM: ' num2str(ppm)]);
    if ( abs(ppm_store(1) - ppm_store(2))/abs(ppm_store(1)) ) > (1/20)
        idx_in_fo_search_set = hit_pss_fo_set_idx(valid_idx);
        
        fo_idx = idx_in_fo_search_set - (pss_idx-1).*length(fo_search_set); % remove duplicated fo
        if fo_idx(1) == fo_idx(2)
            disp(['Discard duplicated frequency idx ' num2str(idx_in_fo_search_set)]);
            idx_in_fo_search_set = idx_in_fo_search_set(1);
        end
        
        f_set = fo_search_set(mod(idx_in_fo_search_set-1, length(fo_search_set)) + 1);
        disp(['Period PPM ' num2str(ppm) 'PPM; f_set ' num2str(f_set./1e3) 'kHz']);
        disp(['Final PSS idx ' num2str(pss_idx)]);
        return;
    end
elseif var(ppm_store) > 0.01
    mean_ppm = mean(ppm_store);
    tmp = abs(ppm_store - mean_ppm);
    [~, max_idx] = max(tmp);
    drop_idx = find(ppm_store==ppm_store(max_idx));
    if length(drop_idx) >= (length(ppm_store)*3/8)
        disp('Too many PPM drops. Will not do it.');
        extra_frequency_flag = 1;
    else
        disp(['Drop PPM: ' num2str(ppm_store(drop_idx))]);
        ppm_store(drop_idx) = [];
        valid_idx(drop_idx) = [];
        ppm_idx = ppm_idx - length(drop_idx);
    end
    ppm = mean(ppm_store);
    disp(['Total ' num2str(ppm_idx) ' freq. idx for PPM: ' num2str(valid_idx)]);
    pss_idx = floor( (hit_pss_fo_set_idx(valid_idx)-1)./length(fo_search_set) ) + 1;
    disp(['Total ' num2str(ppm_idx) '  pss. idx for PPM: ' num2str(pss_idx)]);
    disp(['Average PPM: ' num2str(ppm)]);
end

sum_corr_val = zeros(1, ppm_idx);
for i=1:ppm_idx
    col_idx = valid_idx(i);
    sum_corr_val(i) = sum(hit_corr_val(time_location_invalid_record(:,col_idx)==0,col_idx));
end

[~, max_idx] = max(sum_corr_val);
disp(['Freq. idx for f_set: ' num2str(valid_idx(max_idx))]);

idx_in_fo_search_set = hit_pss_fo_set_idx(valid_idx(max_idx));

if extra_frequency_flag == 1
    extra_valid_idx = valid_idx;
    tmp = zeros(1, ppm_idx);
    tmp(drop_idx) = 1;
    reserve_idx = find(tmp==0);
    if sum(drop_idx==max_idx)
        extra_drop_idx = drop_idx;
    elseif sum(reserve_idx==max_idx)
        extra_drop_idx = reserve_idx;
    else
        disp('Abnormal!');
        return;
    end
    extra_valid_idx(extra_drop_idx) = [];
    extra_ppm_idx = ppm_idx - length(extra_drop_idx);
    
    sum_corr_val = zeros(1, extra_ppm_idx);
    for i=1:extra_ppm_idx
        col_idx = extra_valid_idx(i);
        sum_corr_val(i) = sum(hit_corr_val(time_location_invalid_record(:,col_idx)==0,col_idx));
    end

    [~, max_idx] = max(sum_corr_val);
    disp(['Extra Freq. idx for f_set: ' num2str(extra_valid_idx(max_idx))]);
    pss_idx = floor( (hit_pss_fo_set_idx(extra_valid_idx(max_idx))-1)./length(fo_search_set) ) + 1;
    disp(['Extra  Pss. idx for f_set: ' num2str(pss_idx)]);
    
    idx_in_fo_search_set = [idx_in_fo_search_set hit_pss_fo_set_idx(extra_valid_idx(max_idx))];
end

% add fo from other PSSs if there are
pss_idx = floor( (idx_in_fo_search_set-1)./length(fo_search_set) ) + 1;

extra_pss_set = zeros(1, num_pss);
for idx=1:num_pss
    if ~isempty(find(pss_idx==idx, 1))
        extra_pss_set(idx) = 1;
    end
end
if sum(extra_pss_set)<3
    extra_pss_set = find(extra_pss_set==0);
    tmp_hit_pss_fo_set_idx = hit_pss_fo_set_idx(valid_idx_backup);
    exist_pss_idx = floor( (tmp_hit_pss_fo_set_idx-1)./length(fo_search_set) ) + 1;
    for idx=1:length(extra_pss_set)
        extra_pss_idx = extra_pss_set(idx);
        col_set = find(exist_pss_idx==extra_pss_idx);
        
        if ~isempty(col_set)
            sum_corr_val = zeros(1, length(col_set));
            for i=1:length(col_set)
                col_idx = valid_idx_backup( col_set(i) );
                sum_corr_val(i) = sum(hit_corr_val(time_location_invalid_record(:,col_idx)==0,col_idx));
            end

            [~, max_idx] = max(sum_corr_val);
            disp(['Extra Freq. idx for f_set (multi-PSS): ' num2str(valid_idx_backup( col_set(max_idx) ))]);
            pss_idx = floor( (hit_pss_fo_set_idx(valid_idx_backup( col_set(max_idx) ))-1)./length(fo_search_set) ) + 1;
            disp(['Extra  Pss. idx for f_set (multi-PSS): ' num2str(pss_idx)]);

            idx_in_fo_search_set = [idx_in_fo_search_set hit_pss_fo_set_idx(valid_idx_backup( col_set(max_idx) ))];
        end
    end
end

pss_idx = floor( (idx_in_fo_search_set-1)./length(fo_search_set) ) + 1;
fo_idx = idx_in_fo_search_set - (pss_idx-1).*length(fo_search_set); % remove duplicated fo
[fo_idx, sort_idx] = sort(fo_idx);
idx_in_fo_search_set = idx_in_fo_search_set(sort_idx);
drop_idx = find(diff(fo_idx)==0) + 1;
if ~isempty(drop_idx)
    tmp_idx = zeros(1, 2*length(drop_idx));
    tmp_idx(2:2:end) = drop_idx;
    tmp_idx(1:2:end) = drop_idx-1;
    disp(['Discard duplicated frequency idx (multi-PSS) ' num2str( idx_in_fo_search_set(tmp_idx) )]);
    idx_in_fo_search_set(drop_idx) = [];
end

f_set = fo_search_set(mod(idx_in_fo_search_set-1, length(fo_search_set)) + 1);

disp(['Period PPM ' num2str(ppm) 'PPM; f_set ' num2str(f_set./1e3) 'kHz']);
disp(['Final PSS idx ' num2str(pss_idx)]);

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
