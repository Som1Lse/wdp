#!/bin/sh

echo "k-value,mutation-rate,trials,initialization,fitness,success-rate,sim-time,mutations,mean-degree,mean-p" > results/variable-p-comparisons.csv

# base
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 25 -P datasets/node_probabilities.txt | tee -a results/variable-p-comparisons.csv &

# degree-init, p-fitness
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 25 -P datasets/node_probabilities.txt degree | tee -a results/variable-p-comparisons.csv &

# p-init, degree-fitness
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 25 -P datasets/node_probabilities.txt probability | tee -a results/variable-p-comparisons.csv &

# base
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 50 -P datasets/node_probabilities.txt | tee -a results/variable-p-comparisons.csv &

# degree-init, p-fitness
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 50 -P datasets/node_probabilities.txt degree | tee -a results/variable-p-comparisons.csv &

# p-init, degree-fitness
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 50 -P datasets/node_probabilities.txt probability | tee -a results/variable-p-comparisons.csv &

wait

# base
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 100 -P datasets/node_probabilities.txt | tee -a results/variable-p-comparisons.csv &

# degree-init, p-fitness
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 100 -P datasets/node_probabilities.txt degree | tee -a results/variable-p-comparisons.csv &

# p-init, degree-fitness
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 100 -P datasets/node_probabilities.txt probability | tee -a results/variable-p-comparisons.csv &

# base
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 200 -P datasets/node_probabilities.txt | tee -a results/variable-p-comparisons.csv &

# degree-init, p-fitness
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 200 -P datasets/node_probabilities.txt degree | tee -a results/variable-p-comparisons.csv &

# p-init, degree-fitness
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 200 -P datasets/node_probabilities.txt probability | tee -a results/variable-p-comparisons.csv &

wait

# base
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 250 -P datasets/node_probabilities.txt | tee -a results/variable-p-comparisons.csv &

# degree-init, p-fitness
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 250 -P datasets/node_probabilities.txt degree | tee -a results/variable-p-comparisons.csv &

# p-init, degree-fitness
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 250 -P datasets/node_probabilities.txt probability | tee -a results/variable-p-comparisons.csv &

# base
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 300 -P datasets/node_probabilities.txt | tee -a results/variable-p-comparisons.csv &

# degree-init, p-fitness
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 300 -P datasets/node_probabilities.txt degree | tee -a results/variable-p-comparisons.csv &

# p-init, degree-fitness
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 300 -P datasets/node_probabilities.txt probability | tee -a results/variable-p-comparisons.csv &

wait

# base
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 400 -P datasets/node_probabilities.txt | tee -a results/variable-p-comparisons.csv &

# degree-init, p-fitness
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 400 -P datasets/node_probabilities.txt degree | tee -a results/variable-p-comparisons.csv &

# p-init, degree-fitness
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 400 -P datasets/node_probabilities.txt probability | tee -a results/variable-p-comparisons.csv &

# top-degree
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 25 -P datasets/node_probabilities.txt -d | tee -a results/variable-p-comparisons.csv &

# bot-p
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 25 -P datasets/node_probabilities.txt -e | tee -a results/variable-p-comparisons.csv &

# top-degree
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 50 -P datasets/node_probabilities.txt -d | tee -a results/variable-p-comparisons.csv &

# bot-p
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 50 -P datasets/node_probabilities.txt -e | tee -a results/variable-p-comparisons.csv &

wait

# top-degree
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 100 -P datasets/node_probabilities.txt -d | tee -a results/variable-p-comparisons.csv &

# bot-p
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 100 -P datasets/node_probabilities.txt -e | tee -a results/variable-p-comparisons.csv &

# top-degree
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 200 -P datasets/node_probabilities.txt -d | tee -a results/variable-p-comparisons.csv &

# bot-p
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 200 -P datasets/node_probabilities.txt -e | tee -a results/variable-p-comparisons.csv &

wait

# top-degree
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 250 -P datasets/node_probabilities.txt -d | tee -a results/variable-p-comparisons.csv &

# bot-p
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 250 -P datasets/node_probabilities.txt -e | tee -a results/variable-p-comparisons.csv &

# top-degree
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 300 -P datasets/node_probabilities.txt -d | tee -a results/variable-p-comparisons.csv &

# bot-p
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 300 -P datasets/node_probabilities.txt -e | tee -a results/variable-p-comparisons.csv &

# top-degree
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 400 -P datasets/node_probabilities.txt -d | tee -a results/variable-p-comparisons.csv &

# bot-p
build/moran Explore-Update/datasets/vk/vk.txt -s 500 -k 400 -P datasets/node_probabilities.txt -e | tee -a results/variable-p-comparisons.csv &

wait
