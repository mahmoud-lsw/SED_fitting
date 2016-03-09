!#/bin/bash

python fitspec_2comp_zpt_UDS.py 0 255 &
python fitspec_2comp_zpt_UDS.py 256 510 &
python fitspec_2comp_zpt_UDS.py 511 765 &
python fitspec_2comp_zpt_UDS.py 766 1019 &

wait

python zpt_merge_2comp_UDS.py
python output_merge_2comp_UDS.py

wait

#rm flux_ratios_2comp_0_467.txt flux_ratios_2comp_468_935.txt flux_ratios_2comp_936_1402.txt flux_ratios_2comp_1403_1870.txt
#rm photoz_2comp_0_467.txt photoz_2comp_468_935.txt photoz_2comp_936_1402.txt photoz_2comp_1403_1870.txt
