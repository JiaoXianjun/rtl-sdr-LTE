% Jiao Xianjun (putaoshu@msn.com; putaoshu@gmail.com)
% Find out PSS location.
% A script of project: https://github.com/JiaoXianjun/rtl-sdr-LTE

function [hit_pss_fo_set_idx, hit_time_idx, hit_corr_val] = pss_moving_corr(s, pss_fo_set, fo_search_set)

len_pss = size(pss_fo_set, 1);
num_fo_pss = size(pss_fo_set, 2);
hit_pss_fo_set_idx = 0;
hit_time_idx = 0;
hit_corr_val = 0;

len = length(s);

corr_store = zeros((len - (len_pss-1)), num_fo_pss);

for i=1:(len - (len_pss-1))
    chn_tmp = s(i:(i+len_pss-1));
    tmp = abs(chn_tmp'*pss_fo_set).^2;
    corr_store(i,:) = tmp;
end

corr_store_shifts = cell(1,3);
pss_period = [(19200/2)-1, (19200/2), (19200/2)+1];
num_half_radioframe = floor( (len - (len_pss-1))/pss_period(2) );
max_peak_all = zeros(1, 3*num_fo_pss);
max_idx_all = zeros(1, 3*num_fo_pss);
for i=1:3
    corr_store_tmp = corr_store(1:pss_period(i), : );
    for j=2:num_half_radioframe
        sp = (j-1)*pss_period(i) + 1;
        ep = j*pss_period(i);
        corr_store_tmp = corr_store_tmp + corr_store(sp:ep, : );
    end
    corr_store_tmp = corr_store_tmp + circshift(corr_store_tmp, [-1,0]) + circshift(corr_store_tmp, [1,0]);
    corr_store_shifts{i} = corr_store_tmp;
    sp = (i-1)*num_fo_pss + 1;
    ep = i*num_fo_pss;
    [max_peak_all(sp:ep), max_idx_all(sp:ep)] = max(corr_store_tmp, [], 1);
end

[max_peak, sort_idx] = sort(max_peak_all, 'descend');
max_reserve = 8;
ppm = zeros(1, max_reserve);
f_set = zeros(1, max_reserve);
shift_set = [-1 0 1];
for i=1:max_reserve
    shift_idx = floor( ( sort_idx(i)-1 )/num_fo_pss ) + 1; % 1 -- -1; 2 -- 0; 3 -- 1
    fo_pss_idx = sort_idx(i) - (shift_idx-1)*num_fo_pss;
    
    % calculate frequency offset
    f_set(i) = fo_search_set(mod(fo_pss_idx-1, length(fo_search_set)) + 1);
    
    % calculate PPM
    combined_seq = corr_store_shifts{shift_idx};
    combined_seq = combined_seq(:, fo_pss_idx);
    peak_val = max_peak(i);
end

% % % % -----------------------------old ---------------------------
% len_pss = size(pss_fo_set, 1);
% num_fo_pss = size(pss_fo_set, 2);
% hit_pss_fo_set_idx = 0;
% hit_time_idx = 0;
% hit_corr_val = 0;
% 
% len = length(s);
% % metric_record = zeros(len - (len_pss-1), num_fo_pss);
% 
% len_half_store = 64;
% corr_store = zeros(2*len_half_store+1, num_fo_pss);
% end_idx = inf;
% for i=1:(len - (len_pss-1))
%     chn_tmp = s(i:(i+len_pss-1));
%     
% %     chn_tmp = chn_tmp - mean(chn_tmp);
%     
%     chn_tmp = sqrt(len_pss).*chn_tmp./sqrt( sum(abs(chn_tmp).^2) ); % normalize
%     
%     tmp = abs(chn_tmp'*pss_fo_set).^2;
% %     metric_record(i,:) = tmp;
% 
%     corr_store(2:end,:) = corr_store(1:(end-1),:);
%     corr_store(1,:) = tmp;
%     
%     if sum(tmp>th)
%         current_idx = i;
%         end_idx = current_idx + len_half_store;
%         break;
%     end
% end
% 
% % num_fo = num_fo_pss/3;
% % figure;
% % for i=1:20
% %     subplot(4,5,i); plot(metric_record(:,i));
% % end
% % figure;
% % for i=1:20
% %     subplot(4,5,i); plot(metric_record(:,1*num_fo+i));
% % end
% % figure;
% % for i=1:20
% %     subplot(4,5,i); plot(metric_record(:,2*num_fo+i));
% % end
% % 
% % figure;
% % for i=21:40
% %     subplot(4,5,i-20); plot(metric_record(:,i));
% % end
% % figure;
% % for i=21:40
% %     subplot(4,5,i-20); plot(metric_record(:,1*num_fo+i));
% % end
% % figure;
% % for i=21:40
% %     subplot(4,5,i-20); plot(metric_record(:,2*num_fo+i));
% % end
% % % return;
% 
% 
% if end_idx ~= inf
%     last_idx = min(end_idx, (len - (len_pss-1)));
% 
%     for i=(current_idx+1):last_idx
%         chn_tmp = s(i:(i+len_pss-1));
% 
%         chn_tmp = sqrt(len_pss).*chn_tmp./sqrt( sum(abs(chn_tmp).^2) ); % normalize
% 
%         tmp = abs(chn_tmp'*pss_fo_set).^2;
%     %     metric_record(i,:) = tmp;
% 
%         corr_store(2:end,:) = corr_store(1:(end-1),:);
%         corr_store(1,:) = tmp;
%     end
%     [max_val, max_idx] = max(corr_store, [], 1);
%     [max_val, sort_idx] = sort(max_val, 'descend');
%     
%     tmp = max_val - (max_val(1)/2);
%     idx = find(tmp<0, 1, 'first');
%     if isempty(idx)
%         num_valid = num_fo_pss;
%     else
%         num_valid = idx-1;
%     end
%     hit_pss_fo_set_idx = sort_idx(1:num_valid);
%     hit_corr_val = max_val(1:num_valid);
%     hit_time_idx = last_idx - max_idx(hit_pss_fo_set_idx) + 1;
% end
% 
