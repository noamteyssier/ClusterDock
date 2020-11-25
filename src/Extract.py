#!/usr/bin/env python3

import argparse
import glob
import os
import re
import sys
import gzip

import pandas as pd
from multiprocessing import Pool
from tqdm import tqdm

class Directory:
    def __init__(self, dirname):
        self.dirname = dirname
        self.subclusters = glob.glob(
            os.path.join(self.dirname, "subcluster*")
        )

        self.re_header = re.compile("^[ ]+mol#")
        self.re_values = re.compile("^[ ]+[0-9]")
        self.re_ignore = re.compile(
            "^ 9[0-9]+ 9[0-9]+|skip_size|colors|poses|<|>|no_match|bump|clashes"
            )

        self.re_first_ws = re.compile("^[ ]+")
        self.re_ws = re.compile("[ ]+")

        self.groups = []

        self.raw_mol2_data = []

    def process_line(self, line):
        line = self.re_first_ws.sub("", line)
        line = self.re_ws.sub("\t", line)
        tokens = line.strip().split("\t")

        tokens[0] = int(tokens[0])

        for i in range(2, len(tokens)):
            if '.' in tokens[i]:
                tokens[i] = float(tokens[i])
            else:
                tokens[i] = int(tokens[i])

        return tokens

    def process_outdock(self, fn, sc_id):
        f = open(fn, "r")

        cluster_id=0

        while True:
            try:
                line = next(f)

                if self.re_header.search(line):
                    cluster_id += 1

                elif self.re_ignore.search(line):
                    continue

                elif self.re_values.search(line):
                    tokens = [sc_id, cluster_id] + self.process_line(line)
                    self.groups.append(tokens)

            except StopIteration:
                break

    def process_subcluster(self, sc):
        sc_id = int(sc[-4:])
        outdock_fn = os.path.join(sc, "OUTDOCK")
        self.process_outdock(outdock_fn, sc_id)

    def mol2_parse_hit(self, line, f):
        lines = [line] + [next(f) for _ in range(28)]
        return lines

    def read_mol2(self, sub_idx, cls_idx):
        fn = os.path.join(
            self.dirname,
            "subcluster{:04}".format(sub_idx),
            "test.{:04}.mol2.gz".format(cls_idx)
        )

        f = gzip.open(fn, "rt")

        re_lookup = [
            re.compile("###[ ]+Name:[ ]+{}".format(i))
            for i in self.mol_lookup[(sub_idx, cls_idx)]
        ]

        found = 0

        while True:
            try:
                line = next(f)
                for r in re_lookup:
                    if r.search(line):
                        data = self.mol2_parse_hit(line, f)
                        self.raw_mol2_data.append([sub_idx, cls_idx] + data)
                        found += 1

            except StopIteration:
                break

    def parse_mol2(self):
        def prep_coords(ox_set):
            return [[i for i in ox.strip().split(" ") if i != ""][-3:] for ox in ox_set]

        frame = []
        for d in self.raw_mol2_data:

            sub_idx = d[0]
            cls_idx = d[1]
            name = d[2].split(" ")[-1].strip()
            oxr_set = prep_coords(d[12:16])
            oxs_set = prep_coords(d[16:20])
            total = round(float(d[-1].split(" ")[-1].strip()), 2)

            for i in range(4):
                row = [sub_idx, cls_idx, name, float(total)] + \
                    oxr_set[i] + oxs_set[i]
                frame.append(row)


        self.mol2_frame = pd.DataFrame(
            frame,
            columns = \
                ["sub_idx", "cls_idx", "mol_name", "Total"] + \
                ["rec_{}".format(x) for x in ["x", "y", "z"]] + \
                ["lig_{}".format(x) for x in ["x", "y", "z"]]
            )

    def process_mol2(self):
        to_process = self.uniq_scores.\
            loc[:,['sub_idx', 'cls_idx', 'mol_name']].\
            drop_duplicates().\
            values

        self.mol_lookup = {}

        for i in to_process:
            sub_idx, cls_idx, mol_name = i
            file = (sub_idx, cls_idx)
            if file not in self.mol_lookup:
                self.mol_lookup[file] = []

            self.mol_lookup[file].append(mol_name)

        for f in self.mol_lookup:
            sub_idx, cls_idx = f
            self.read_mol2(sub_idx, cls_idx)



        self.parse_mol2()

        # select for hits
        self.mol2_hits = self.uniq_scores.\
            merge(
                self.mol2_frame,
                left_on = ['sub_idx', 'cls_idx', 'Total', "mol_name"],
                right_on = ['sub_idx', 'cls_idx', 'Total', "mol_name"]
                ).\
            loc[:, self.mol2_frame.columns]

    def write(self):

        # for enrichment
        self.scores.drop('cls_idx', axis = 1).to_csv(
            os.path.join(self.dirname, "extract_all.sort.txt"),
            sep = "\t", index=False, header=None
        )

        # for enrichment
        self.uniq_scores.drop('cls_idx', axis = 1).to_csv(
            os.path.join(self.dirname, "extract_all.sort.uniq.txt"),
            sep = "\t", index=False, header=None
        )


        # for downstream processing
        self.uniq_scores.to_csv(
            os.path.join(self.dirname, "unique_scores.tab"),
            sep = "\t", index=False
        )

        # for mol2 coordinates
        self.mol2_hits.to_csv(
            os.path.join(self.dirname, "unique_coords.tab"),
            sep = "\t", index=False
        )


    def run(self):

        # collect OUTDOCK information
        [self.process_subcluster(s) for s in self.subclusters]

        # aggregate scores into dataframe
        self.scores = pd.DataFrame(
            self.groups,
            columns = [
                "sub_idx", "cls_idx", "mol_idx", "mol_name", "flex_code",
                "matched", "nscored", "time", "hac", "setnum", "matnum",
                "rank", "cloud", "elect", "gist", "vdW", "psol", "asol",
                "inter", "rec_e", "rec_d", "r_hyd", "Total"
                ]
            ).sort_values("Total")
        self.scores = self.scores[~self.scores.mol_name.isna()]

        # take best scoring of unique molecules
        self.uniq_scores = self.scores.\
            sort_values("Total").\
            drop_duplicates("mol_name")


        # collect mol2 data
        self.process_mol2()


        # write data
        self.write()

def run_dir(d):
    dir = Directory(d)
    dir.run()

def get_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "-i", '--input', required=True, nargs = "+",
        help = "input directory(s) to process"
        )
    p.add_argument(
        "-n", "--threads", required=False,
        default = 8, type=int,
        help = 'number of threads to process on'
    )
    args = p.parse_args()

    return args

def main():
    args = get_args()

    p = Pool(args.threads)
    p.map(run_dir, tqdm(args.input, desc='Processing Directories'))
    p.close()


if __name__ == '__main__':
    main()
