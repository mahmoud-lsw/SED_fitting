python model_mags_make_2c.py

wait

cd zeropoints/
./fitspec_2comp_zpt.sh

wait

python MC_offsets.py

wait

cd ..
./fitspec_2comp.sh
