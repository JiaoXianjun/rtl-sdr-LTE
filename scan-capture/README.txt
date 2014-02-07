frequency-2575-2595MHz and frequency-2635-2655MHz contain captured I&Q data at target center frequency.

After they are decompressed, you will see:

*.it are captured by CellSearch (https://github.com/Evrytania/LTE-Cell-Scanner)
*.bin are captured by rtl_sdr (http://sdr.osmocom.org/trac/wiki/rtl-sdr)

Three it files (0000, 0001, 0002) are from low exact and high frequencies.

Example rules of bin file name:

f647_s1.92_g0_10s.bin

    frequency -- 647+1998MHz
sampling rate -- 1.92MHz
         gain -- 0 (automatic gain)
          10s -- length of capture period

