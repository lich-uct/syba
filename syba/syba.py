#!/usr/bin/env python

"""
Library for classification of chemical structures based on frequency of fragments in training sets.
Originaly developed for predicting synthetic feasibility, but possibly can be used for different class problem.
"""

# Copyright (c) 2019 Milan Vorsilak
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import math
import os
import gzip
from rdkit.Chem import AllChem as Chem


class SybaClassifier:
    """
    SYBA clasifier
    """
    def __init__(self, neighbourhood=4):
        self.fragments = {}
        self.ALL_OFF_FRAGS_SCORE = None
        self.pNS = None
        self.neighbourhood = neighbourhood
        self.chiral_scores = {
                 0: (1.8084159515032894, -0.4675609913440742),
                 1: (0.8694169068939861, -0.25951468034597774),
                 2: (-0.10727513159632096, 0.025893146969645207),
                 3: (-1.4655318562675832, 0.19407129808044743),
                 4: (-2.939230216688627, 0.1883123534101861),
                 5: (-4.536353950205819, 0.2060140915608036)
                }


    def fitDefaultScore(self):
        this_dir, this_filename = os.path.split(__file__)
        with gzip.open(os.path.join(this_dir, "resources", "syba4.csv.gz"), mode="rt") as counts:
            self.fitFromCountFile(counts)


    def fitFromCountFile(self, reader):
        """
        """
        header = reader.readline()  # header is skipped
        spls = reader.readline().strip().split(",")
        n_syn_structures = int(spls[1])
        n_non_structures = int(spls[2])
        self.pNS = math.log((n_non_structures + 2) / (n_syn_structures + 2))
        self.ALL_OFF_FRAGS_SCORE = 0.0
        self.fragments = {}
        for line in reader:
            spls = line.strip().split(",")
            frg = int(spls[0])
            counts = (int(spls[1]), int(spls[2]))
            scores = calculateScore(counts[0], counts[1], n_syn_structures, n_non_structures)
            self.fragments[frg] = scores
            self.ALL_OFF_FRAGS_SCORE+=self.pNS + scores[1]


    def fitFromScoreFile(self, reader, pNS=0):
        """
        """
        header = reader.readline()  # header is skipped
        spls = reader.readline().strip().split(",")
        self.pNS = pNS
        self.ALL_OFF_FRAGS_SCORE = 0.0
        self.fragments = {}
        for line in reader:
            spls = line.strip().split(",")
            frg = int(spls[0])
            scores = (float(spls[1]), float(spls[2]))
            self.fragments[frg] = scores
            self.ALL_OFF_FRAGS_SCORE+=self.pNS + scores[1]


    def predict(self, smi=None, mol=None, useChiral=True):
        if smi:
            mol = Chem.MolFromSmiles(smi)
        frgs = Chem.GetMorganFingerprint(mol,self.neighbourhood)
        score = self.ALL_OFF_FRAGS_SCORE
        for frg,y in frgs.GetNonzeroElements().items():
            if frg in self.fragments:
                frg_score = self.fragments[frg]
                score -= frg_score[1]
                score += frg_score[0]
        if useChiral:
            score += self._getChiralScore(mol, includeAbsentChiralities=True)
        return score


    def predict_with_only_present_fragments(self, smi=None, mol=None, useChiral=True):
        if smi:
            mol = Chem.MolFromSmiles(smi)
        frgs = ch.GetMorganFingerprint(mol,self.neighbourhood)
        score = 0.0
        for frg,y in frgs.GetNonzeroElements().items():
            if x in self.fragments:
                score += self.pNS + self.fragments[frg][0]
        if useChiral:
            score += self._getChiralScore(mol, includeAbsentChiralities=False)
        return score


    def fragment_contribution(self, smi=None, mol=None):
        if smi:
            mol = Chem.MolFromSmiles(smi)
        contribs = [0 for a in mol.GetAtoms()]
        info={}
        fp = Chem.GetMorganFingerprint(mol,self.neighbourhood,bitInfo=info)
        for key, vals in info.items():
            if key in self.fragments:
                s = self.fragments[key][0]
                for a,n in vals:
                    contribs[a]+=s
        return contribs

    def _getChiralScore(self, mol, includeAbsentChiralities=True):
        chiral_count = Chem.CalcNumAtomStereoCenters(mol)
        if chiral_count > max(self.chiral_scores.keys()):
            chiral_count = max(self.chiral_scores.keys())
        elif chiral_count < 0:
            raise Exception("Compound can't have negative chiral count")
        if includeAbsentChiralities:
            return sum(c[0] if k == chiral_count else c[1] for k, c in self.chiral_scores.items())
        else:
            return self.chiral_scores[chiral_count][0]


def SmiMolSupplier(reader, header=False, smi_col=0, delim=","):
    if header:
        header = reader.readline().strip().split(delim)
    for i,line in enumerate(reader):
        spls = line.strip().split(delim)
        smi = spls[smi_col]
        try:
            yield (Chem.MolFromSmiles(smi), *spls)
        except Exception as e:
            print("SMILES at line {} is not readable. '{}'".format(i+1, line.strip()))
            raise e


def getEnvironment(m, values):
    for value in values:
        atom = value[0]
        radius = value[1]
        env = Chem.FindAtomEnvironmentOfRadiusN(m,radius+1,atom)
        amap={}
        submol=Chem.PathToSubmol(m,env,atomMap=amap)
        if len(amap) > 0:
            try:
                return Chem.MolToSmiles(submol,rootedAtAtom=amap[atom])
            except KeyError as e:
                print(m.GetAtomWithIdx(atom).GetSymbol())
                print(amap)
                raise e
        else:
            # print(amap, m.GetAtomWithIdx(atom).GetSymbol(), Chem.MolToSmiles(m))
            return m.GetAtomWithIdx(atom).GetSymbol()


def getFragmentSmiles(m, neighbourhood=2):
    info={}
    fp = Chem.GetMorganFingerprint(m,neighbourhood,bitInfo=info)
    env = []
    for i in info:
        env.append((i,getEnvironment(m,info[i])))
    return env


def processFile(mol_supplier, get_fragments=False, neighbourhood=4):
    """
    For given compressed file, function returns number of all compounds and all Morgan fragment types with corresponding number of compounds
    in which they are found.

    path to data set (gzip csv file)
    smi_col is column where are SMILES, first column has index 0!
    """
    fs = {}
    if get_fragments:
        def updateFragments(m):
            for frag, smi in getFragmentSmiles(m):
                d = fs.setdefault(frag,[0, set()])
                d[0] += 1
                d[1].add(smi)
    else:
        def updateFragments(m):
            for frag in Chem.GetMorganFingerprint(m,neighbourhood).GetNonzeroElements().keys():
                d = fs.setdefault(frag, [0])
                d[0] += 1
    n_of_compounds = 0
    for m in mol_supplier:
        updateFragments(m)
        n_of_compounds += 1
    return fs, n_of_compounds


def mergeFragmentCounts(feasible_fragments, infeasible_fragments):
    fragments = set(feasible_fragments.keys())
    fragments.update(infeasible_fragments.keys())
    fragment_counts = {}

    for f in fragments:
        fragment_counts[f] = (feasible_fragments.get(f, [0])[0], infeasible_fragments.get(f, [0])[0])
    return fragment_counts


def calculateScore(feasible_count, infeasible_count, feasible_compound_count, infeasible_compound_count):
    return math.log((feasible_count+1)/(infeasible_count+1)), math.log((feasible_compound_count-feasible_count+1)/(infeasible_compound_count-infeasible_count+1))


def writeScoreFile(filename, fragments_counts, compound_counts):
    with open(filename, "w") as out:
        out.write("fragment_id,score_presence,score_absence\n")
        for fragment_id, counts in fragments_counts.items():
            out.write("{},{},{}\n".format(fragment_id, *calculateScore(counts[0], counts[1], compound_counts[0], compound_counts[1])))


def writeCountFile(filename, fragments_counts, compound_counts):
    with open(filename, "w") as out:
        out.write("fragment_id,feasible_compounds,infeasible_compounds\n")
        out.write("total_compounds,{},{}\n".format(compound_counts[0], compound_counts[1]))
        for fragment_id, counts in fragments_counts.items():
            out.write("{},{},{}\n".format(fragment_id, counts[0], counts[1]))


def main():
    import argparse
    import sys

    parser = argparse.ArgumentParser(description="Syba classifier scores structures, more positive values are for easier- and more negative for harder-to-synthesize.")

    parser.add_argument('input_file', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help="path to an csv, otherwise STDIN is used. Expected file format is 'ID,SMILES' and with header.")
    parser.add_argument('output_file', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help="path to output file, otherwise STDOUT is used.")
    args = parser.parse_args()
    header = True
    
    this_dir, this_filename = os.path.split(__file__)
    syba = SybaClassifier()
    syba.fitDefaultScore()

    with args.output_file as out, args.input_file as inp:
        if header:
            inp.readline()
        supp = SmiMolSupplier(inp, smi_col=1)
        out.write("ID,SMILES,SYBA\n")
        for mol, *spls in supp:
            out.write(f"{spls[0]},{Chem.MolToSmiles(mol)},{syba.predict(mol=mol)}\n")


if __name__ == "__main__":
    main()
