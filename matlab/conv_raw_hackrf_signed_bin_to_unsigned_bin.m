function conv_raw_hackrf_signed_bin_to_unsigned_bin(bin_filename)

fid = fopen(bin_filename, 'r');

if fid==-1
    disp('conv_raw_hackrf_signed_bin_to_unsigned_bin: Can not open file for read!');
    return;
end

tmp_store = fread(fid, inf, 'int8');
fclose(fid);

fid = fopen(bin_filename, 'w');

if fid==-1
    disp('conv_raw_hackrf_signed_bin_to_unsigned_bin: Can not open file for write!');
    return;
end

tmp_store = tmp_store + 128;
fwrite(fid, tmp_store, 'uint8');

fclose(fid);

