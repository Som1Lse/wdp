#!/bin/sh

echo "accuracy,probability,MCtime,mutations" > results/required-simulations.csv
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 100 -p 1.0 -a 100 -i | tee -a results/required-simulations.csv
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 100 -p 1.0 -a 1000 -i | tee -a results/required-simulations.csv
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 100 -p 1.0 -a 10000 -i | tee -a results/required-simulations.csv
