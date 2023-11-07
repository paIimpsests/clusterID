"""
   plot_hist.py
------------------
Script to plot histogram from output of 'clusterID' script (default is a "logs.txt" file).
Default is the raw count (not normalized) with integer bins, this can easily be tweaked by looking at 'numpy.histogram' documentation.
Can save the histogram with desired format if '-o <output>' option is specified.

OPTIONS:
 -i [str]: SOURCE file
 -o [str]: OUTPUT file
 -h: print help

"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import getopt


def parse_input(argv):
    arg_input = "logs.txt"
    arg_output = ""
    arg_help = "{0} -i <input> -o <output> -h".format(argv[0])

    try:
        opts, args = getopt.getopt(argv[1:], "hi:o:", ["help", "input=", "output="])
    except:
        print(arg_help)
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(arg_help)  # print the help message
            sys.exit(2)
        elif opt in ("-i", "--input"):
            arg_input = arg
        elif opt in ("-o", "--output"):
            arg_output = arg

    return arg_input, arg_output


SOURCE, OUTPUT = parse_input(sys.argv)

# Fetch data from SOURCE
data = []
with open(SOURCE, "r") as f:
    line_count = 0
    next_header = 1
    for line in f:
        line_count += 1
        if line_count == next_header:
            next_header += int(line.split("&")[1]) + 1
        else:
            size, count = np.array(line.split("\t"), dtype=np.float32)
            for _ in range(int(count)):
                data.append(size)
f.close()
data = np.array(data, dtype=np.int32)

# Plot histogram
fig, ax = plt.subplots()
n, bins, patches = ax.hist(
    data,
    bins=np.max(data) - 1,
    stacked=False,
    density=False,
    align="left",
    rwidth=0.9,
)
ax.set(
    xlabel="cluster size",
    ylabel="count",
    xlim=[0.5, None],
)

# Save histogram if OUTPUT is specified
if OUTPUT != "":
    plt.savefig(OUTPUT, format=OUTPUT.split(".")[-1])

plt.show()
