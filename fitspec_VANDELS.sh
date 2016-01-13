!#/bin/bash

python fitspec_VANDELS.py 0 25 &
python fitspec_VANDELS.py 26 50 &
python fitspec_VANDELS.py 51 75 &
python fitspec_VANDELS.py 76 100 &

wait

python output_merge_2comp_V.py

wait

#rm photoz_V_0_13.txt photoz_V_14_25.txt photoz_V_26_37.txt photoz_V_38_50.txt

#wait

#python plotspec.py
