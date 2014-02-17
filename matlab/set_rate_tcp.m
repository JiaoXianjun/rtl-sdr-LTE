% Jiao Xianjun (putaoshu@msn.com; putaoshu@gmail.com)
% parameter setting via tcp to rtl_tcp
% A script of project: https://github.com/JiaoXianjun/multi-rtl-sdr-calibration

function tcp_obj = set_rate_tcp(tcp_obj, rate)
fwrite(tcp_obj, 2, 'uint8');
fwrite(tcp_obj, uint32(rate), 'uint32');
