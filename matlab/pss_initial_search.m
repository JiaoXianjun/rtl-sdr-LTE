% Jiao Xianjun (putaoshu@msn.com; putaoshu@gmail.com)
% Find out PSS location.
% A script of project: https://github.com/JiaoXianjun/rtl-sdr-LTE

function [hit_flag, hit_idx, hit_sss_avg_power, hit_pss_power, pss_idx] = pss_initial_search(s, pss, th)
num_pss = size(pss, 2);
len_pss = size(pss, 1);

hit_flag = false;
hit_idx = -1;
pss_idx = -1;
hit_sss_avg_power = inf;
hit_pss_power = inf;

mv_len = 128+9; % the last OFDM symbol (SSS) of slot 1 and 11
distance_pss_sss = 9+10+9 + 3*128; % normal CP case
store_for_moving_avg = 999.*ones(num_pss, distance_pss_sss);
sampling_rate = 1.92e6; % LTE spec
fo_search_set = [-100e3:5e3:100e3];
num_fo = length(fo_search_set);

pss_set = zeros(len_pss, num_fo, num_pss);
for i=1:num_pss
    pss_set(:,:,i) = kron(ones(1, num_fo), pss(:,i)).*exp(1i.*2.*pi.*(1./sampling_rate).*(0:(len_pss-1)).'*fo_search_set);
end

sum_power = sum(store_for_moving_avg(:, (end-mv_len+1):end), 2);

len = length(s);
metric_record = zeros(len - (len_pss-1), num_fo, num_pss);

% calculate matrix for segment corr
len_seg = len_pss;
% len_seg = 32;
num_seg = floor(len_pss/len_seg);
len_corr = len_seg*num_seg;
fetch_matrix = zeros(num_seg, len_corr);
for i=1:num_seg
    sp = (i-1)*len_seg + 1;
    ep = sp + len_seg - 1;
    fetch_matrix(i, sp:ep) = 1;
end

current_power = zeros(num_pss, 1);
for i=1:(len - (len_pss-1))
    chn_tmp = s(i:(i+len_pss-1));
    
    chn_tmp = sqrt(len_pss).*chn_tmp./sqrt( sum(abs(chn_tmp).^2) ); % normalize
    
    for j=1:num_pss
        pss = pss_set(:,:,j);
        tmp = conj(pss).*kron(ones(1, num_fo), chn_tmp);
        tmp = fetch_matrix*tmp(1:len_corr,:);
        tmp = sum( abs(tmp).^2, 1 );
        metric_record(i,:,j) = tmp;
        current_power(j) = max(tmp);
%         metric_record(i,1,j) = current_power(j);
    end
    
    current_power_dB = 10.*log10(current_power);
    avg_power_dB = 10.*log10(sum_power./mv_len);
    peak_to_avg = current_power_dB - avg_power_dB;

    logic_tmp = (peak_to_avg > th);
    if sum(logic_tmp)
        hit_flag = true;
        pss_idx = find(logic_tmp==1);
        hit_idx = i;
        hit_pss_power = current_power_dB(pss_idx);
        hit_sss_avg_power = hit_pss_power - peak_to_avg(pss_idx);
        break;
    else
        sum_power = sum_power - store_for_moving_avg(:, end);
        sum_power = sum_power + store_for_moving_avg(:, end-mv_len);
        
        store_for_moving_avg(:, 2:end) = store_for_moving_avg(:, 1:(end-1));
        store_for_moving_avg(:, 1) = current_power;
    end
    
end

% figure;
% subplot(3,1,1); plot(metric_record(:,1,1));
% subplot(3,1,2); plot(metric_record(:,1,2));
% subplot(3,1,3); plot(metric_record(:,1,3));

figure;
for i=1:20
    subplot(4,5,i); plot(metric_record(:,i,1));
end
figure;
for i=1:20
    subplot(4,5,i); plot(metric_record(:,i,2));
end
figure;
for i=1:20
    subplot(4,5,i); plot(metric_record(:,i,3));
end

figure;
for i=21:40
    subplot(4,5,i-20); plot(metric_record(:,i,1));
end
figure;
for i=21:40
    subplot(4,5,i-20); plot(metric_record(:,i,2));
end
figure;
for i=21:40
    subplot(4,5,i-20); plot(metric_record(:,i,3));
end

