% Jiao Xianjun (putaoshu@msn.com; putaoshu@gmail.com)
% parameter setting via tcp to rtl_tcp
% A script of project: https://github.com/JiaoXianjun/multi-rtl-sdr-calibration


function tcp_obj = set_gain_tcp(tcp_obj, gain)

if gain
    fwrite(tcp_obj, 3, 'uint8');
    fwrite(tcp_obj, uint32(1), 'uint32');
    fwrite(tcp_obj, 4, 'uint8');
    fwrite(tcp_obj, uint32(gain), 'uint32');
else
    fwrite(tcp_obj, 3, 'uint8');
    fwrite(tcp_obj, uint32(0), 'uint32');
end

