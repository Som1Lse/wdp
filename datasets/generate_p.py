# code to generate independent, random P values for each node
import os, random

src = "./"
max = 0.7
min = 0.1
size = 500

def gen_P():
    return random.uniform(min,max)


def IO():
    with open(os.path.join(src, "node_probabilities.txt"), 'w') as fout:
        for x in range(0, size-1):
            fout.write('{:f}'.format(gen_P()) + "\n")

if __name__ == "__main__":
    IO()
