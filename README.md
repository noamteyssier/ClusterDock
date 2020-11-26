# Matching Sphere Cluster Analysis

These are scripts for performing matching sphere clustering and for analyzing the downstream results.


# Prerequisites

Clone the ClusterDock repository and initialize scripts. This will make sure everything is pointing to the right place.

```bash
# I recommend putting this in a place you can access easily (~/bin, ~/lib)
git clone https://github.com/noamteyssier/dock_cluster

cd dock_cluster
./init.sh
```

# Preparing Input Directory
The initial clustering script requires a particular format to run correctly. 
It requires you have all the necessary files already prepared (see blastermaster tutorial)

It requires a directory as an argument : this is the meta directory of the receptor ready to be docked.

```
# Input Directory : 
meta/
	dockfiles/
		ligand.desolv.heavy
		ligand.desolv.hydrogen
		matching_spheres.sph
		trim.electrostatics.phi
		vdw.bmp
		bwd.parms.amb.mindock
		vdw.vdw
	INDOCK
	ligands.names
	decoys.names
	enrichment_sdi
```

# Preparing Clusters
The script to prepare clusters is found under :
`<git_path>/src/ClusterSpheres.py`

The script will run k-means on the matching sphere set and generate a directory for each k for as many iterations as it is given. The match_goal will be scaled to (Match_Goal * 1/k) if the flag is provided to the run. The num_sdi flag controls how many SDI clusters each run will be split into - the current default is 20 as is the default with blastermaster

It will make a dirlist including all directories and subdirectories for easy submission on wynton + gimel

```
usage: ClusterSpheres.py [-h] -i INPUT -k NUM_CLUSTERS [NUM_CLUSTERS ...] [-n NUM_ITER] [-m] [-f] [-v] [-s NUM_SDI]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input meta directory to prepare for clustering
  -k NUM_CLUSTERS [NUM_CLUSTERS ...], --num_clusters NUM_CLUSTERS [NUM_CLUSTERS ...]
                        Number of clusters to split matching sphere set into (multiple arguments allowed)
  -n NUM_ITER, --num_iter NUM_ITER
                        Number of iterations for each cluster k to create
  -m, --scale_match_goal
                        scale match goal in INDOCK files to reflect 1/k matches
  -f, --overwrite       overwrite created directories
  -v, --verbose         increase verbosity
  -s NUM_SDI, --num_sdi NUM_SDI
                        Number of clusters to split sdi set into

```

Here are some example usages : 
```bash
# prepare a clustered run with 7 subclusters
./ClusterSpheres.py -i meta/ -k 7 -n 1

# prepare a clustered run of 5 subclusters with 20 different k-means runs
./ClusterSpheres.py -i meta/ -k 5 -n 20

# prepare clustered runs for k2 and k4 with 15 different k-means runs each
./ClusterSpheres.py -i meta/ -k 2,4 -n 15

# prepare clustered runs for k1 through k5 with 10 different k-means runs each
./ClusterSpheres.py -i meta/ -k {2..5} -n 10

# prepare clustered runs for k1 through k5 with 10 different k-means runs each, overwriting old ones
./ClusterSpheres.py -i meta/ -k {2..5} -n 10 -f

# prepare an unclustered run
./ClusterSpheres.py -i meta/ -k 1 -n 1

# prepare clustered runs for k3 with 10 different k-means runs with scaled match_goal parameter
./ClusterSpheres.py -i meta/ -k 3 -n 10 -m

```

# Running DOCK on each cluster

The current release of DOCK doesn't include the matching sphere usage output for each molecule so you will need to run the branch I am running. This is included in the repository alongside its submission wrapper scripts. 

`<git_path>/bin/dock64`


```bash
# to submit an individual DOCK run
<git_path>/src/dock64

# to submit a DOCK run with queue
<git_path>/src/submit.csh

```

# Extracting Results

DOCK results are a bit annoying to parse out, and I wanted to be flexible to existing scripts out there but also write code that makes dataframes that are easy to work with downstream. I have written a script that extracts all relevant data for cluster analysis - I decided to write it in Julia instead of python because it was roughly the same development time with a 10-100x speedup in terms of data processing - especially once parallel processing was involved. The only cost is that there is a long startup time to this script as it reads in the precompiled packages. You will need a version of Julia >1 for this script. I recommend preparing all the packages beforehand. The script is found here `<git_path>/src/Extract.jl`. This tutorial will assume you have julia on your path.

## Prerequisite : Prepare Julia

```bash
# check version (version at time of writing : 1.4.1)
julia --version

# enter REPL
julia
```

Once you're in the REPL just drop this code in. 
It'll take a little while to run if it's your first time building these. 

```julia
# This is a good time to make a coffee 

using Pkg
Pkg.add(["ArgParse", "Distributed", "DataFrames", "Statistics", "Printf", "GZip", "CSV"])
using ArgParse, Distributed, DataFrames, Statistics, Printf, GZip, CSV
```
## Output Files

This script operates on each cluster directory it is given via the `-i` flag which accepts as many arguments as is possible from the command line. 

*You will have poor performance launching it independently for each one*

Each directory will have 4 files written :
1) extract_all.sort.txt 
2) extract_all.sort.uniq.txt
3) time_and_enrichment.tab
4) coords.tab

(1) and (2) are written more for legacy purposes and to be backwards compatible with other scripts that may use these - they will not be used again in downstream analyses.

(3) is the number of time spent in each SDI subcluster alongside the AUC and LogAUC of the run. This is calculated using the best scoring pose of each molecule.

(4) is the coordinates of the molecules considered hits. Essentially they are the coordinates of the molecules in extract_all.sort.uniq.txt that are below a certain threshold of Total energy. This threshold is determined given a quantile (the -q flag) but the default is 0.1. This will only accept molecules as hits if their energy is in the bottom 10% of all unique molecule energy. 

## Example Usages
```bash

# run on a single clustered directory (k2_1)
julia Extract.jl -i k2_1

# run on all k2 directories
julia Extract.jl -i k2*/

# run on all k directories and use 16 threads
julia -p 16 Extract.jl -i k*/

# run on all k directories for all receptors using 24 threads
julia -p 24 Extract.jl -i {AA2AR,EGFR,AMPC}/k*/ 

# run on all k directories using 12 threads with a quantile of 0.2
julia -p 12 Extract.jl -i k*/ -q 0.2
```

