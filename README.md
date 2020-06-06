This is the code from our (I, [loke202](https://github.com/loke202) and [cookiepig](https://github.com/cookiepig)) bachelor's thesis.

To compile, run
```sh
mkdir build
cd build
cmake ..
make
```

Aside from a C++14 compiler and standard library, the code has no dependencies.

To run make, sure you clone with `--recursive` to get the VK dataset from https://github.com/nd7141/Explore-Update.

The shell scripts can be used to produce the files in the `results` folder.

`datasets/node_probabilities.txt` was generated using the `generate_p.py` in the same directory. Unfortunately it contains a typo, where it only generates 499 entries instead of 500, which causes the code to be incorrect. We believe the results of this are relatively minor due to how the p-values were used (primarily in comparisons, so they don't propagate throught the code). Unfortunately, there is not enough time to fix this, and we felt it would be better to document everything exactly as is, instead of hiding the mistake.
