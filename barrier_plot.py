import numpy as np
import matplotlib.pyplot as plt

filename = "DGn.txt"
n, DGn = np.loadtxt(filename, unpack=True, delimiter="\t")
fig, ax = plt.subplots()
ax.scatter(n, DGn, marker=".")
ax.set(xlabel=r"$\textrm{cluster size }n$",
       ylabel=r"$\beta DG(n)$",
       title=r"$\textrm{{{title}}}$".format(title=filename),
       )
plt.show()

