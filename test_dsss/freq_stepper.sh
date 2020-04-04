#!/bin/bash

# im folgendem soll per hackrf-one ein Frequenzbereich auf moegliche 
# Spreizsignale geprueft werden 

# "hackrf_transfer" is von osmocom
# test_dsss it hausgemacht

# Antenneneingang nicht unbeschalten lassen

freq_start=1000000
freq_stop=6000000
freq_step=50000


for((freq=$freq_start; $freq<$freq_stop; freq=$freq+$freq_step))
do
    echo "check freq: ${freq}"
    
    # "-r -" stdout, "-f" Mittenfrequenz", "-s" Samplerate, "-n" Samples
    # "-l", if gain
    hackrf_transfer -r - -f ${freq} -s 1000000 -n 65536 -l 24 | ./test_dsss
done

echo "accomplished"
