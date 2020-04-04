#!/bin/bash

# im folgendem soll per hackrf-one ein Frequenzbereich auf moegliche 
# Spreizsignale geprueft werden 

# "hackrf_transfer" is von osmocom
# test_dsss it hausgemacht

freq_start=1250000
freq_stop=599750000
freq_step=500000


for((freq=$freq_start; $freq<$freq_stop; freq=$freq+$freq_step))
do
    echo "check freq: ${freq}"
    
    hackrf_transfer -r /dev/null -f ${freq} -n 65536 -l 24 | ./test_dess
done

echo "accomplished"
