% Jiao Xianjun (putaoshu@msn.com; putaoshu@gmail.com)
% convert a sequence to a mat. each column shift previos to left by 1
% for example input sequence 1 2 3 4 5 6 7, and len_body=4
% output matrix:
% 1 2 3 4
% 2 3 4 5
% 3 4 5 6
% 4 5 6 7

% A script of project: https://github.com/JiaoXianjun/multi-rtl-sdr-calibration

function r = lin2col_shift_mat(s, len_body)

if length(s) < len_body
    disp('Length of input is too short!');
    r = -1;
    return;
end

len_tail = length(s) - len_body;

r = toeplitz(s, [s(1) zeros(1, len_tail)]);
r = r((len_tail+1):end, end:-1:1);
