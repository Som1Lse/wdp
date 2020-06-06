#!/bin/sh
# to warm your room at night try
# echo "k-value,mutation-rate,trials,initialization,fitness,success-rate,sim-time,mutations,mean-degree,mean-p" > results/edge-weight-comparisons.csv; echo {1..25} | xargs -n 1 -P your_cores ./run_moran_edge_weights_tests.sh

# base
echo "$1"
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k "$1" -E  | tee -a results/edge-weight-comparisons.csv &

# top-degree
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k "$1" -E -d | tee -a results/edge-weight-comparisons.csv &
echo "$1"
