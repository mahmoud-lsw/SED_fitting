!#/bin/bash

python fitspec.py 0 &
python fitspec.py 3 &
python fitspec.py 6 &
python fitspec.py 9 &

wait

python output_merge.py

