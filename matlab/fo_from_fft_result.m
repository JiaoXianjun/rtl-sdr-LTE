% Jiao Xianjun (putaoshu@msn.com; putaoshu@gmail.com)
% estimation frequency offset from fft power spectrum
% A script of project: https://github.com/JiaoXianjun/rtl-sdr-LTE

function [fo, int_phase_rotate] = fo_from_fft_result(spectrum_in, sampling_rate)
% each column of spectrum_in is an power spectrum sequence

fft_len = size(spectrum_in,1);
spectrum_in = [spectrum_in( ((fft_len/2)+1):end, : ); spectrum_in( 1:(fft_len/2), : )];
[~, max_idx] = max(spectrum_in, [], 1);
int_phase_rotate = 2.*pi.*(max_idx - ((fft_len/2) + 1 ) )./fft_len;
fo = sampling_rate.*int_phase_rotate./(2*pi);
