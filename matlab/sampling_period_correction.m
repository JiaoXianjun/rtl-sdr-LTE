% Jiao Xianjun (putaoshu@msn.com; putaoshu@gmail.com)

% Convert input signal(with sampling period shorter than expectation by
% ppm_val PPM) to output signal(sampling period error corrected)

% A script of project: https://github.com/JiaoXianjun/rtl-sdr-LTE

function r = sampling_period_correction(s, ppm_val)
sampling_time_in = ( 0 : (length(s)-1) ).*(1-ppm_val*1e-6);
sampling_time_out = 0 : (length(s)-1);
if ppm_val > 0
    ep =  floor( (length(s)-1).*(1-ppm_val*1e-6) );
    sampling_time_out = sampling_time_out(1: ep);
end
r = interp1(sampling_time_in, s, sampling_time_out, 'cubic').';
