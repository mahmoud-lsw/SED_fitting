!#/bin/bash

python fitspec_2comp.py 0 13 &
python fitspec_2comp.py 14 25 &
python fitspec_2comp.py 26 37 &
python fitspec_2comp.py 38 50 &

wait

python output_merge_2comp.py

wait

#rm photoz_0_13.txt photoz_14_25.txt photoz_26_37.txt photoz_38_50.txt
