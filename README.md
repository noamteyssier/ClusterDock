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

The current release of DOCK doesn't include the matching sphere usage output for each molecule
so you will need to run the branch I am running. This is included in the repository alongside
its submission wrapper scripts. 

`<git_path>/bin/dock64`


```bash
# to submit an individual DOCK run
<git_path>/src/dock64

# to submit a DOCK run with queue
<git_path>/src/submit.csh

```


