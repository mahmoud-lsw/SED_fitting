!#/bin/bash

python fitspec_2comp_zpt.py 0 467 &
python fitspec_2comp_zpt.py 468 935 &
python fitspec_2comp_zpt.py 936 1402 &
python fitspec_2comp_zpt.py 1403 1870 &

wait

python zpt_merge_2comp.py
python output_merge_2comp.py

wait

rm flux_fracoffsets_2comp_0_467.txt flux_fracoffsets_2comp_468_935.txt flux_fracoffsets_2comp_936_1402.txt flux_fracoffsets_2comp_1403_1870.txt
rm photoz_2comp_0_467.txt photoz_2comp_468_935.txt photoz_2comp_936_1402.txt photoz_2comp_1403_1870.txt
