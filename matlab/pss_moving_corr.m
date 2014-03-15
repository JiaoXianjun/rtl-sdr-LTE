% Jiao Xianjun (putaoshu@msn.com; putaoshu@gmail.com)
% Find out PSS location.
% A script of project: https://github.com/JiaoXianjun/rtl-sdr-LTE

function [hit_pss_fo_set_idx, hit_time_idx, hit_corr_val] = pss_moving_corr(s, pss_fo_set, th)

len_pss = size(pss_fo_set, 1);
num_fo_pss = size(pss_fo_set, 2);
hit_pss_fo_set_idx = 0;
hit_time_idx = 0;
hit_corr_val = 0;

len = length(s);
% metric_record = zeros(len - (len_pss-1), num_fo_pss);

len_half_store = 64;
corr_store = zeros(2*len_half_store+1, num_fo_pss);
end_idx = inf;
for i=1:(len - (len_pss-1))
    chn_tmp = s(i:(i+len_pss-1));
    
%     chn_tmp = chn_tmp - mean(chn_tmp);
    
    chn_tmp = sqrt(len_pss).*chn_tmp./sqrt( sum(abs(chn_tmp).^2) ); % normalize
    
    tmp = abs(chn_tmp'*pss_fo_set).^2;
%     metric_record(i,:) = tmp;

    corr_store(2:end,:) = corr_store(1:(end-1),:);
    corr_store(1,:) = tmp;
    
    if sum(tmp>th)
        current_idx = i;
        end_idx = current_idx + len_half_store;
        break;
    end
end

% num_fo = num_fo_pss/3;
% figure;
% for i=1:20
%     subplot(4,5,i); plot(metric_record(:,i));
% end
% figure;
% for i=1:20
%     subplot(4,5,i); plot(metric_record(:,1*num_fo+i));
% end
% figure;
% for i=1:20
%     subplot(4,5,i); plot(metric_record(:,2*num_fo+i));
% end
% 
% figure;
% for i=21:40
%     subplot(4,5,i-20); plot(metric_record(:,i));
% end
% figure;
% for i=21:40
%     subplot(4,5,i-20); plot(metric_record(:,1*num_fo+i));
% end
% figure;
% for i=21:40
%     subplot(4,5,i-20); plot(metric_record(:,2*num_fo+i));
% end
% % return;


if end_idx ~= inf
    last_idx = min(end_idx, (len - (len_pss-1)));

    for i=(current_idx+1):last_idx
        chn_tmp = s(i:(i+len_pss-1));

        chn_tmp = sqrt(len_pss).*chn_tmp./sqrt( sum(abs(chn_tmp).^2) ); % normalize

        tmp = abs(chn_tmp'*pss_fo_set).^2;
    %     metric_record(i,:) = tmp;

        corr_store(2:end,:) = corr_store(1:(end-1),:);
        corr_store(1,:) = tmp;
    end
    [max_val, max_idx] = max(corr_store, [], 1);
    [max_val, sort_idx] = sort(max_val, 'descend');
    
    tmp = max_val - (max_val(1)/2);
    idx = find(tmp<0, 1, 'first');
    if isempty(idx)
        num_valid = num_fo_pss;
    else
        num_valid = idx-1;
    end
    hit_pss_fo_set_idx = sort_idx(1:num_valid);
    hit_corr_val = max_val(1:num_valid);
    hit_time_idx = last_idx - max_idx(hit_pss_fo_set_idx) + 1;
end

