#!/bin/sh

echo "mutation-percent,probability,MCtime,mutations" > results/mutation-rates.csv
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -p 1.0 -m 0.05 | tee -a results/mutation-rates.csv
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -p 1.0 -m 0.1 | tee -a results/mutation-rates.csv
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -p 1.0 -m 0.2 | tee -a results/mutation-rates.csv
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -p 1.0 -m 0.3 | tee -a results/mutation-rates.csv
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -p 1.0 -m 0.4 | tee -a results/mutation-rates.csv
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -p 1.0 -m 0.5 | tee -a results/mutation-rates.csv
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -p 1.0 -m 0.6 | tee -a results/mutation-rates.csv
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -p 1.0 -m 0.7 | tee -a results/mutation-rates.csv
