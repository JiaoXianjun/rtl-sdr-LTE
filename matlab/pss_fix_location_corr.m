% Jiao Xianjun (putaoshu@msn.com; putaoshu@gmail.com)
% Perform pss correlation at specific locations
% A script of project: https://github.com/JiaoXianjun/rtl-sdr-LTE

function [hit_time_idx, max_val] = pss_fix_location_corr(s, i_set, pss_fo_set, hit_pss_fo_set_idx)
len_pss = size(pss_fo_set, 1);

sp = i_set(1);
ep = i_set(2) + len_pss - 1;
r = lin2col_shift_mat(s(sp:ep), len_pss);

% r = r - kron( ones(len_pss,1), mean(r, 1) );

r = sqrt(len_pss).*r./kron( ones(len_pss,1), sqrt( sum(abs(r).^2, 1) ) ); % normalize

pss_fo_specific_set = pss_fo_set(:, hit_pss_fo_set_idx);

corr_val = abs(r'*pss_fo_specific_set).^2;

[max_val, max_idx] = max(corr_val, [], 1);

hit_time_idx = sp + max_idx - 1;
