#!/usr/bin/env python3

import numpy as np
import pandas as pd
import argparse
import sys
import os

from multiprocessing import Pool

class SelectHits:
    def __init__(self, dirname, percentile=10, include_ligands=False, label=None):
        self.dirname = dirname
        self.percentile = percentile
        self.include_ligands = include_ligands
        self.label = label

        self.coord_fn = os.path.join(self.dirname, "unique_coords.tab")
        self.output_fn = os.path.join(self.dirname, "top_coords.tab")

        self.validate_dir()
        self.run()

    def validate_dir(self):
        if not os.path.exists(self.coord_fn):
            sys.exit()

    def select_top(self):
        pc = np.percentile(self.coord_frame.Total, [self.percentile])[0]
        self.coord_frame = self.coord_frame[self.coord_frame.Total < pc]

    def exclude_ligands(self):
        self.coord_frame = self.coord_frame.\
            drop(columns = ['lig_x', 'lig_y', 'lig_z'])

    def run(self):
        self.coord_frame = pd.read_csv(self.coord_fn, sep="\t")

        self.select_top()
        if not self.include_ligands:
            self.exclude_ligands()

        if self.label:
            self.coord_frame['label'] = self.label

        self.coord_frame.to_csv(
            self.output_fn, sep="\t", index=False
        )

def get_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "-i", "--input", nargs='+',
        required=True, type = str,
        help = "Input directory to process"
    )
    p.add_argument(
        "-p", "--percentile", default=10,
        required=False, type = int,
        help = "Percentile to consider a hit (default = 10)"
    )
    p.add_argument(
        "-v", "--include_ligands", action='store_true',
        help = "include ligand coordinates (default = exclude)"
    )
    p.add_argument(
        '-l', '--label',
        required=False, type = str,
        help = 'Label to include on output dataframe'
    )
    p.add_argument(
        '-n', '--num_threads', default=8,
        required=False, type = int,
        help = 'number of threads to pool'
    )
    p.add_argument(
        '-f', '--file',
        required=False, type = str,
        help = "read arguments from file (will ignore -i/--input flag)"
    )
    args = p.parse_args()

    return args

def main():
    args = get_args()

    if args.file:
        dirlist = []
        with open(args.file, "r") as f:
            for line in f:
                dirlist.append(line.strip())

        arglist = [
            [i, args.percentile, args.include_ligands, i] for i in dirlist
        ]

        p = Pool(args.num_threads)
        p.starmap(SelectHits, arglist)
        p.close()

    elif len(args.input) == 1:
        sh = SelectHits(
            args.input[0], args.percentile,
            args.include_ligands, args.label
            )
            
    else:
        arglist = [
            [i, args.percentile, args.include_ligands, i] for i in args.input
        ]
        p = Pool(args.num_threads)
        p.starmap(SelectHits, arglist)
        p.close()



if __name__ == '__main__':
    main()
