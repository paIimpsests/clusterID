"""
   plot_barrier.py
---------------------
Script to plot free energy barrier from output of 'clusterID' script (default is a "DGn.txt" file).
Can save the plot with desired format if '-o <output>' option is specified.

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
    arg_input = "DGn.txt"
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


# Fetch data from SOURCE
SOURCE, OUTPUT = parse_input(sys.argv)
n, DGn = np.loadtxt(SOURCE, unpack=True, delimiter="\t")

# Plot
fig, ax = plt.subplots()
ax.scatter(n, DGn, marker=".")
ax.set(
    xlabel="cluster size n$",
    ylabel="beta DG(n)",
)

# Save plot if OUTPUT is specified
if OUTPUT != "":
    plt.savefig(OUTPUT, format=OUTPUT.split(".")[-1])

plt.show()
