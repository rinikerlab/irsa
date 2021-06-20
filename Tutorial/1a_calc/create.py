import numpy as np
import matplotlib.pyplot as py
import rdkit
from rdkit import *
from rdkit.Chem import *
from rdkit.Chem.rdDistGeom import EmbedMultipleConfs
from rdkit.Chem.rdmolfiles import *
f=open('SMILES','r')
string=str(f.readline())
smile=string
solute_mol = AddHs(MolFromSmiles(smile))
EmbedMultipleConfs(solute_mol, numConfs=250, clearConfs=True, pruneRmsThresh=0.5, numThreads=8)
for i in range(250):
    rdkit.Chem.rdmolfiles.MolToXYZFile(solute_mol, str(i)+".xyz", confId=i)


