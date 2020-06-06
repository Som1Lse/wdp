#!/bin/bash

Run(){
    build/im                                                                                  \
        Explore-Update/datasets/vk/vk.txt Explore-Update/datasets/vk/vk_mv.txt -b 10          \
        Explore-Update/datasets/vk/vk_mem.txt Explore-Update/datasets/vk/vk_seeds.txt -s 250  \
        simpath -t 0.00001 sigmoid_degree "$@"
}

Pids=()

Run top -k 60 nodes                >results/mv-s250-top-nodes.csv &

Run top -k 60 edges                >results/mv-s250-top-edges.csv &

Run top -k 60 simpath              >results/mv-s250-top-simpath.csv &

Run greedy -k 60 mc -i 100         >results/mv-s250-greedy-mc-1.csv &

Run greedy -k 60 mc -i 1000        >results/mv-s250-greedy-mc-2.csv &

# This takes eight hours+ to run. Feel free to uncomment it if you don't need to use your computer for a while.
# Run greedy -k 60 mc -i 10000       >results/mv-s250-greedy-mc-3.csv &

Run greedy -k 60 simpath -t 0.1    >results/mv-s250-greedy-simpath-1.csv &

Run greedy -k 60 simpath -t 0.01   >results/mv-s250-greedy-simpath-2.csv &

Run greedy -k 60 simpath -t 0.001  >results/mv-s250-greedy-simpath-3.csv &

Run greedy -k 60 simpath -t 0.0001 >results/mv-s250-greedy-simpath-4.csv &

Run greedy -k 60 simpath -t 0.5    >results/mv-s250-greedy-simpath-5.csv &

wait
