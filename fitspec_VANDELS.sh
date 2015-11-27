!#/bin/bash

python fitspec_VANDELS.py 0 13 &
python fitspec_VANDELS.py 14 25 &
python fitspec_VANDELS.py 26 37 &
python fitspec_VANDELS.py 38 50 &

wait

python output_merge_2comp_V.py

wait

rm photoz_V_0_13.txt photoz_V_14_25.txt photoz_V_26_37.txt photoz_V_38_50.txt

wait

python plotspec.py
