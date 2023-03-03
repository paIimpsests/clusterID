# clusterID.c

Cluster identification routine based on ten Wolde's local orientational order parameter to differentiate between liquid- and solid-like particles. The routine produces a `.sph` file where clusters are highlighted (see below for details on supported file formats). Can also produce a log of identified clusters as well as free energy calculations for the formation of a cluster of a given size. The latter is only relevant for SOURCE file consisting of a large number of system snapshots with a constant number of particles. These can be obtained from experimental data or brute-force nucleation simulations (Molecular Dynamics, unbiased Monte-Carlo).

Current implementation is still *kind of wacky* for experimental data submitted in `.xyz` format as it relies on a larger-then-expected system size to bypass periodic boundary conditions (PBC).

## Build

`cd` in workind directory and `make` to build.

## Usage

Run `clusterID` with `-h` to get the list of available options (also listed below).

Standard usage is to run `clusterID` with any `[OPTIONS]` followed by the `SOURCE` file which supported formats are listed below.

Some example `SOURCE` files can be found under the `example` folder.

## Options

 * `-b` **[double]**: ten Wolde's local orientational order parameter dot-product cutoff for two particles to be identified as connected
 * `-n` **[int]**: number of connections cutoff for a particle to be identified as solid-like
 * `-o` **[double]**: ten Wolde's local orientationnal order parameter dot-product cutoff for two particles to be identified as part of the same cluster
 * `-r` **[double]**: radius cutoff for nearest neighbors (NN) identification. Automatically triggers choice of cutoff radius method for NN identification
 * `-L`: only highlights the largest identified cluster
 * `-f` **[int]**: choice of SOURCE file type among:
 
   **0** for `.sph` files
        
   **1** for `.xyz` files (default)

 * `-m` **[int]**: choice of NN identification method among:
 			
   **0** for SANN method (default)
      
   **1** for 12NN method
      
   **2** for cutoff radius method
 
 * `-l`: saves cluster logs to a file
 * `-G`: calculates free energy for the formation of a cluster of size n, saves it to a file
 * `-h`: displays help and exit

## Supported SOURCE file formats

Only two file formats are supported at the moment: `.sph` and `.xyz`. A source file can contain multiple snapshots, consecutively typed. The format for each snapshot depends on file format and is presented below. 

### `.sph` snapshot format

Consists of $N+2$ lines where $N$ is the number of particles:

```
&N
L L L
a x y z r
```

And the third line is repeated $N$ times. $L$ is the 1D size of the cubic system box. $a$ (or any letter) indicates particle type. $x$, $y$, and $z$ are the 3D coordinates of the particle. $r$ is the particle radius.

### `.xyz` snapshot format

Consists of $N+2$ lines where $N$ is the number of particles.

```
N
*******
A x y z
```

And the third lines is repeated $N$ times.  $x$, $y$, and $z$ are the 3D coordinates of the particle. It is assumed that they are normalized in terms of particle diameter (hence $r = 0.5$ is assumed).

Note that the current implementation calculated the system box size to be much larger to escape PBC for this supplied file format.

## Analysis details

WIP

## To Do List

- [ ] Allow to chose between boundary conditions (PBC for simulated data vs. hard-wall for experimental data)

