rtl-sdr-LTE
===========

run: matlab/CellSearch.m

have fun!

News:

2014-03-16

All scripts are ported into C/C++ in https://github.com/JiaoXianjun/LTE-Cell-Scanner.
Now you can use the stand alone C/C++ program to scan both TDD and FDD LTE Cell, both without external LNB and with external LNB!

----

Play with LTE signal captured by rtl-sdr.

Introduction video(inside China): http://v.youku.com/v_show/id_XNjc1MjIzMDEy.html

Introduction video(outside China): http://www.youtube.com/watch?v=4zRLgxzn4Pc

By a enthusiast <putaoshu@msn.com> <putaoshu@gmail.com> of Software Defined Radio.

Some files are borrowed from https://github.com/Evrytania/LTE-Cell-Scanner, https://github.com/Evrytania/Matlab-Library and https://github.com/JiaoXianjun/multi-rtl-sdr-calibration

As China has announced TD-LTE deployment officially, I want to decode it with LTE-Cell-Scanner (written by James Peroulas: https://github.com/Evrytania/LTE-Cell-Scanner).

Then there are some discoveries:

Frequencies of LTE in China are discovered by searching in internet:

FDD band:

China Mobile:  1880-1900MHz (maybe it is also TDD)

China telecom: 1755-1785MHZ(uplink)/1850-1880MHz(downlink)

China unicom:  1955-1980MHz(uplink)/2145-2170MHz(downlink)

TDD band:

China Mobile:  2575-2635 MHz

China telecom: 2635-2655 MHz

China unicom:  2555-2575 MHz

Then those bands are scanned with my MATLAB scanner: https://github.com/JiaoXianjun/multi-rtl-sdr-calibration (see README_for_spectrum_scanner.txt in that repo. Or just search rtl-sdr in http://www.mathworks.com/matlabcentral/fileexchange/)

See those .png files for spectrum plots. (For the band above 2.5GHz which has exceeded range of rtl-sdr dongle, a MMDS LNB is used to extend the band to above 2.5GHz.
See MMDS-LNB-LO1998-to-extend-dongle-band.jpg. I learn this method from http://blog.cyberexplorer.me/2014/01/sniffing-and-decoding-nrf24l01-and.html and https://github.com/omriiluz/NRF24-BTLE-Decoder.
The LO of my LNB is 1998MHz. It means that when the dongle is tuned to 600MHz, it actually receives 600+1998Mhz!)

LTE-Cell-Scanner decodes LTE MIB successfully in 1850-1880MHz band, but unsuccessful for other bands even they seems pretty like LTE 20MHz spectrum.

I think it maybe because relative time location of PSS and SSS is different in TDD from FDD. It is confirmed with James Peroulas.

Inspired by James Peroulas, initial exploration with only PSS correlation is done on captured IQ data, and valid PSS has been seen there! (Try it with matlab/test_td_lte_pss.m)

Now TD-LTE signal is identified successfully, maybe adding TD-LTE support to LTE-Cell-Scanner is a good idea (Done! See: https://github.com/JiaoXianjun/LTE-Cell-Scanner).

But now LTE-Cell-Scanner only works for situation where there isn't MMDS LNB, because LTE-Cell-Scanner assumes analytic relationship between carrier and sampling frequency.
When there is external LNB, the relationship won't be maintained anymore. I offer a work around. A module sampling_ppm_f_search_set_by_pss.m is used to do
sampling frequency estimation (get k_factor). Thus LTE-Cell-Scanner can work only for uncertain carrier frequency (By getting k_factor from sampling_ppm_f_search_set_by_pss). By this way
the whole scanner drops analytic relationship between carrier and sampling frequency. See CellSearch.m for detailed informaiton.

Please join me if you are also interested in this. Please see TODO firstly.

News
=======================
2014-03-16

All scripts are ported into C/C++ in https://github.com/JiaoXianjun/LTE-Cell-Scanner.
Now you can use the stand alone C/C++ program to scan both TDD and FDD LTE Cell, both without external LNB and with external LNB!

2014-01-19:

CellSearch.m works for both TDD and FDD mode under both with LNB and without LNB.

2014-01-07:

PSS has been detected successfully at frequency 2645MHz! See comments in matlab/test_td_lte_pss.m, and then run it.

Usage
=======================
2014-01-19:

Make sure your rtl-sdr dongle works fine (http://sdr.osmocom.org/trac/wiki/rtl-sdr). Then run "rtl_tcp -p 1234 -d 0" in shell. Then:

CellSearch.m (Open it. See comments. And try!)

2014-01-07:

--test_td_lte_pss.m

Enter matlab directory, and run:

test_td_lte_pss('../scan-capture/frequency-2635-2655MHz/f2645_s1.92_g0_1s.bin')

--extract_part_from_rtl_sdr_bin.m

If you have huge captured bin file, extract_part_from_rtl_sdr_bin.m can be used to extract part of it as a new bin file.

Contributing
=======================
You are welcome to send pull requests in order to improve the project.

See TODO list included in the source distribution first (If you want).


