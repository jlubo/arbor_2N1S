#!/bin/bash

# plot averaged data in png format with all existing samples
python3 plotAveragedResults.py -protocol=STET -out_format=png -sampling=1
python3 plotAveragedResults.py -protocol=WTET -out_format=png -sampling=1
python3 plotAveragedResults.py -protocol=SLFS -out_format=png -sampling=1
python3 plotAveragedResults.py -protocol=WLFS -out_format=png -sampling=1


