#!/bin/sh

echo "k-value,probability,MCtime,mutations" > results/test.csv
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 10 -p 1.0 -a 1000 -r 1 -g | tee -a results/test.csv
