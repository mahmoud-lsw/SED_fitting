!#/bin/bash

python fitspec.py 0 20 &
python fitspec.py 21 42 &
python fitspec.py 43 64 &
python fitspec.py 65 126 &

wait

python output_merge.py

wait

#rm photoz_0_20.txt photoz_21_42.txt photoz_43_64.txt photoz_65_126.txt
