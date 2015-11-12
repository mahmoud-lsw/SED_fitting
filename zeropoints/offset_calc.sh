!#/bin/bash

python offset_calc.py 0 467 &
python offset_calc.py 468 935 &
python offset_calc.py 936 1402 &
python offset_calc.py 1403 1870 &

wait

python offset_merge.py

wait

rm flux_fracoffsets_0_467.txt flux_fracoffsets_468_935.txt flux_fracoffsets_936_1402.txt flux_fracoffsets_1403_1870.txt
