% Jiao Xianjun (putaoshu@msn.com; putaoshu@gmail.com)
% Extract part of rtl-sdr captured bin into a new file.
% A script of project: https://github.com/JiaoXianjun/rtl-sdr-LTE

function extract_part_from_rtl_sdr_bin(input_file_name, num_skip_bytes, num_extract_bytes, output_file_name)

fid = fopen(input_file_name);

if fid==-1
    disp('Can not open input file!');
    return;
end

s = fread(fid, inf, 'uint8');

fclose(fid);

if (num_skip_bytes+num_extract_bytes)>length(s)
    disp('num_skip_bytes + num_extract_bytes > len_input');
    return;
end

s = s((num_skip_bytes+1) : (num_skip_bytes+num_extract_bytes));

fid = fopen(output_file_name, 'w');

if fid==-1
    disp('Can not open output file!');
    return;
end

count = fwrite(fid, s, 'uint8');

fclose(fid);

if count ~= length(s)
    disp('Number of written bytes is not as expectation!');
end
