% Jiao Xianjun (putaoshu@msn.com; putaoshu@gmail.com)
% Find out best timing period in frequency-time grids
% A script of project: https://github.com/JiaoXianjun/rtl-sdr-LTE

function [first_valid_idx, last_valid_idx] = find_best_timing_location(diff_time_location, len_th)
first_valid_idx = -1;
last_valid_idx = -1;

[len_time, num_f] = size(diff_time_location);
diff_time_location = diff_time_location>=4;
sp = zeros(1, num_f);
ep = zeros(1, num_f);

for i=1:num_f
    tmp = diff_time_location(:,i);
    for len = len_time:-1:len_th
        for j=1:(len_time-len+1)
            a = sum(tmp(j:(j+len-1)));
            if a == 0
                break;
            end
        end
        if a == 0
            sp(i) = j;
            ep(i) = sp(i) + len - 1;
            break;
        end
    end
end

if sum(sp) == 0
    return;
end

first_valid_idx = sp;
last_valid_idx = ep+1;

% max_sp = max(sp);
% min_ep = min(ep(ep>0));
% 
% if min_ep-max_sp <9
%     return;
% end
% 
% first_valid_idx = max_sp;
% last_valid_idx = min_ep+1;
