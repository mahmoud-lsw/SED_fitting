!#/bin/bash

python fitspec_2comp.py 0 1 &
python fitspec_2comp.py 2 3 &
python fitspec_2comp.py 4 5 &
python fitspec_2comp.py 6 7 &

wait

python output_merge_2comp.py

wait

rm photoz_0_1.txt photoz_2_3.txt photoz_4_5.txt photoz_6_7.txt
