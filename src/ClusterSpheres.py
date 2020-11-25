#!/usr/bin/env python3

import numpy as np
import pandas as pd
import argparse
import glob
import sys
import os
import shutil
import re
from tqdm import tqdm
from multiprocess import Pool

from sklearn.cluster import KMeans

class SplitSubDir:
    """
    A recreation of setup_db2_zinc15_file_number:
        splits an enrichment_sdi into N evenly split directories
    """

    def __init__(self, sdi_fn, n=20):
        self.sdi_fn = sdi_fn
        self.n = n

        self.split_db = {
            i : [] for i in range(n)
        }

    def prepare(self):
        with open(self.sdi_fn, 'r') as f:
            for idx, line in enumerate(f):
                self.split_db[idx%self.n].append(line)

    def write(self, fn, idx):
        with open(fn, "w+") as f:
            for l in self.split_db[idx]:
                f.write(l)

class ClusterSPH:
    def __init__(self, input_fn, output_fn, k, seed, verbose=False, multi_fn=True):
        self.input_fn = input_fn
        self.output_fn = output_fn
        self.k = k
        self.seed = seed
        self.verbose = verbose
        self.multi_fn = multi_fn

        self.sph_header = []
        self.sph_lines = []
        self.sph_frame = pd.DataFrame()
        self.sph_matrix = np.array([])
        self.cluster_size = np.array([])

    def read_sph(self):
        with open(self.input_fn, "r+") as f:
            while True:
                try:
                    line = next(f)
                    if line[0] != ' ':
                        self.sph_header.append(line)
                    else:
                        self.sph_lines.append(line)
                except StopIteration:
                    self.build_frame()
                    break

    def write_sph(self):
        if self.multi_fn:
            for idx in np.arange(self.k):
                self.write_subcluster(idx)
        else:
            self.write_single_fn()

    def prepare_final_header(self, idx, replace_idx=False):
        line = self.sph_header[-1]
        line = line.split(' ')
        line[-1] = "{}\n".format(self.cluster_size[idx])
        line = " ".join(line)

        if replace_idx:
            line = line.replace("1", str(idx + 1), 1)
        return line

    def write_subcluster(self, idx):
        ofn = self.output_fn + ".{}.sph"
        with open(ofn.format(idx), "w+") as f:
            c_idx = np.where(self.sph_frame.cluster == idx)[0]

            final_header = self.prepare_final_header(idx)
            for line in self.sph_header[:-1]:
                f.write(line)

            f.write(final_header)

            for line in self.sph_lines[c_idx]:
                f.write(line)

    def write_single_fn(self):
        ofn = self.output_fn + ".sph"
        with open(self.output_fn, "w+") as f:

            # write color matching header
            for line in self.sph_header[:-1]:
                f.write(line)

            for k_idx in range(self.k):
                # isolate lines in a cluster
                c_idx = np.where(self.sph_frame.cluster == k_idx)[0]

                # write cluster header
                final_header = self.prepare_final_header(k_idx, replace_idx=True)
                f.write(final_header)

                # write nodes
                for line in self.sph_lines[c_idx]:
                    f.write(line)

    def build_frame(self):
        self.sph_frame = pd.DataFrame([
            [_ for _ in l.strip().split(" ") if _ != ""] for l in self.sph_lines
            ])
        self.matrix = self.sph_frame.iloc[:, [1,2,3]].values
        self.sph_lines = np.array(self.sph_lines)

    def cluster(self):
        km = KMeans(n_clusters = self.k, random_state=self.seed)
        names = km.fit_predict(self.matrix)
        self.sph_frame['cluster'] = names
        self.cluster_size = np.array([
            np.sum(names == i) for i in np.arange(self.k)
            ])

    def log(self):
        print(
            "Number of Clusters : {}".format(self.k)
            )
        for i in np.arange(self.k):
            print(
                "\tspheres in cluster {} : {}".\
                format(i, self.cluster_size[i])
            )

    def run(self):
        self.read_sph()
        self.cluster()
        self.write_sph()

        if self.verbose:
            self.log()

class PrepareClusters:
    def __init__(self, meta_dir, k_list, num_iter=1, scale_match_goal=False, overwrite=False, verbose=False, num_sdi_clusters=20):
        self.meta_dir = meta_dir
        self.k_list = k_list
        self.num_iter = num_iter
        self.scale_match_goal = scale_match_goal
        self.overwrite = overwrite
        self.verbose = verbose
        self.num_sdi_clusters = num_sdi_clusters
        self.pwd = os.getcwd()

        self.ms_fn = os.path.join(meta_dir, "dockfiles/matching_spheres.sph")

        self.to_symlink = []
        self.to_symlink_dockfiles = []

        self.validate_meta()
        if self.verbose:
            self.summarise_input()

        self.ssd = SplitSubDir(
            os.path.join(self.meta_dir, "enrichment_sdi"),
            n = self.num_sdi_clusters
            )
        self.ssd.prepare()

    def summarise_input(self):
        print("Given Meta Directory : {}".format(self.meta_dir))
        print("Given K List : {}".format(self.k_list))
        print("Given Num Iter : {}".format(self.num_iter))

    def validate_meta(self):
        req_fn = ["INDOCK", "ligands.names", "decoys.names", "enrichment_sdi"]

        checks = 0
        for rfn in req_fn:
            f = os.path.join(self.meta_dir, rfn)
            if not os.path.isfile(f):
                print("ERROR : Required File Missing from given meta dir : \n\t{}\n".format(rfn))
                checks += 1
            self.to_symlink.append(f)

        if not os.path.isdir(os.path.join(self.meta_dir, "dockfiles")):
            print("ERROR : dockfiles directory missing from given meta dir")
            checks += 1

        if not os.path.isfile(self.ms_fn):
            print("ERROR : matching spheres missing from given meta dir dockfiles")
            checks += 1

        if checks > 0:
            sys.exit(-1)

        self.to_symlink_dockfiles = glob.glob(
            os.path.join(self.meta_dir, "dockfiles/*")
        )

    def prepare_INDOCK(self, input_fn, output_fn, k):
        lines = []
        with open(input_fn, "r") as f:
            while True:
                try:
                    line = next(f)

                    # replaces match goal if required
                    if ("match_goal" in line) and (self.scale_match_goal):
                        num = int(line.split(' ')[-1])
                        replace_num = int(num / k)
                        line = line.replace(str(num), str(replace_num))

                    # adds line to write list
                    lines.append(line)

                    # appends k_clusters line to INDOCK
                    if "bump_rigid" in line:
                        cluster_line = "".join(
                            ['k_clusters'] + \
                            [" " for _ in range(20)] + \
                            [str(k)] + \
                            ["\n"]
                            )
                        lines.append(cluster_line)

                except StopIteration:
                    break

        with open(output_fn, "w+") as f:
            for l in lines:
                f.write(l)

    def create_directory(self, dir_name):
        if os.path.isdir(dir_name):
            if not self.overwrite:
                sys.exit(
                    "ERROR : Directory exists please delete \n\t{}\n".\
                        format(dir_name)
                    )
            else:
                if self.verbose:
                    print("Overwriting Directory : {}".format(dir_name))
                shutil.rmtree(dir_name)

        os.makedirs(os.path.join(dir_name, "dockfiles"))

    def populate_directory(self, dir_name, k, n):

        dir_dockfiles = os.path.join(dir_name, "dockfiles")

        # symlink meta files to directory
        for fn in self.to_symlink:
            bn = fn.split("/")[-1]
            if bn == "INDOCK":
                self.prepare_INDOCK(
                    fn, os.path.join(dir_name, bn), k
                    )
            else:
                os.symlink(
                    os.path.join(self.pwd, fn),
                    os.path.join(dir_name, bn)
                    )

        # populate dockfiles directory
        for fn in self.to_symlink_dockfiles:
            bn = fn.split("/")[-1]

            if bn == 'matching_spheres.sph':
                continue

            os.symlink(
                os.path.join(self.pwd, fn),
                os.path.join(dir_dockfiles, bn)
            )

        # overwrite matching spheres with clustered set
        # uses n index as seed for random state
        cl = ClusterSPH(
            self.ms_fn, os.path.join(dir_dockfiles, "matching_spheres.sph"),
            k, n, multi_fn=False
        )
        cl.run()

    def prepare_sdi_subclusters(self, dir_name):

        # create sdi subcluster directories
        for i in range(self.num_sdi_clusters):
            subdir = os.path.join(dir_name, "subcluster{:04d}".format(i))
            os.makedirs(subdir)

            # symlink INDOCK
            os.symlink(
                os.path.join(os.path.join(self.pwd, dir_name), "INDOCK"),
                os.path.join(subdir, "INDOCK")
            )

            # write precomputed SDI
            self.ssd.write(
                os.path.join(subdir, "split_database_index"),
                i
            )

        # write subdirectories to dirlist
        with open(os.path.join(dir_name, "dirlist"), "w+") as f:
            for i in range(self.num_sdi_clusters):
                f.write(
                    "./subcluster{:04d}\n".format(i)
                )

    def write_dirlist(self):
        with open("dirlist", "w+") as f:
            for l in glob.glob("k*/subcluster*"):
                f.write("./{}\n".format(l))

    def prepare_directory(self, k, n):
        dir_name = "k{}_{}".format(k, n)
        self.create_directory(dir_name)
        self.populate_directory(dir_name, k, n)
        self.prepare_sdi_subclusters(dir_name)

    def build_clusters(self):
        set_list = [(k, n) for k in self.k_list for n in range(self.num_iter)]

        p = Pool()
        p.starmap(self.prepare_directory, tqdm(set_list))
        p.close()

        self.write_dirlist()
        self.print_face()

    def print_face(self):
        print("")
        print("----------------------------------------")
        print("                  ^^^^^^^^^^             ")
        print("                ^^^^^^^^^^^^^^            ")
        print("               ^^^^^^^^^^^^^^^^            ")
        print("              ^^^^^^^^^^^^^^^^^^            ")
        print("             |     ---     ---  \           ")
        print("             |     <o>   \  <o>  |   YOU'RE ")
        print("            {            \       \          ")
        print("             |         <..>      |    DONE  ")
        print("       0     |                   |            ")
        print("     _(_)    |  \\\\\\\\\\\\\ //////// |          ")
        print("    o __|    |                   |           ")
        print("    o __|    |  |\___________/|  |         ")
        print("    o __|    \   \  \ |_|_| \ / /          ")
        print("    o __|     \   \__\_|_/__// /          ")
        print("   /   /     /  \_     __    _/          ")
        print("  /   /     /     \____\/____/            ")
        print(" /   /     /              /             ")
        print("----------------------------------------")


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "-i", "--input", required=True, type=str,
        help="Input meta directory to prepare for clustering"
    )
    p.add_argument(
        "-k", "--num_clusters", nargs="+", required=True, type=int,
        help="Number of clusters to split matching sphere set into (multiple arguments allowed)"
    )
    p.add_argument(
        "-n", "--num_iter", default=1, required=False, type=int,
        help="Number of iterations for each cluster k to create"
    )
    p.add_argument(
        "-m", "--scale_match_goal", action="store_true", required=False,
        help="scale match goal in INDOCK files to reflect 1/k matches"
    )
    p.add_argument(
        "-f", "--overwrite", action='store_true', required=False,
        help="overwrite created directories"
    )
    p.add_argument(
        "-v", "--verbose", action='store_true', required=False,
        help="increase verbosity"
    )
    p.add_argument(
        "-s", "--num_sdi", default=20, required=False, type=int,
        help="Number of clusters to split sdi set into"
    )
    args = p.parse_args()
    return args

def main():
    args = get_args()

    pc = PrepareClusters(
        meta_dir = args.input,
        k_list = args.num_clusters,
        num_iter = args.num_iter,
        scale_match_goal = args.scale_match_goal,
        overwrite = args.overwrite,
        verbose = args.verbose,
        num_sdi_clusters = args.num_sdi
    )
    pc.build_clusters()

if __name__ == '__main__':
    main()
