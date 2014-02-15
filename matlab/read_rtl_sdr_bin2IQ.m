% Jiao Xianjun (putaoshu@msn.com; putaoshu@gmail.com)
% read specific number of IQ samples from rtl-sdr captured bin file
% A script of project: https://github.com/JiaoXianjun/rtl-sdr-LTE

function s = read_rtl_sdr_bin2IQ(rtl_sdr_bin_filename, num_sample)

fid = fopen(rtl_sdr_bin_filename);

if fid==-1
    disp('read_rtl_sdr_bin2IQ: Can not open input file!');
    return;
end

[s, count] = fread(fid, 2*num_sample, 'uint8');
fclose(fid);

if count < (2*num_sample)
    disp('read_rtl_sdr_bin2IQ: Please ensure number of samples is sufficient in bin file!');
    s = -1;
    return;
end

s = raw2iq(s);
