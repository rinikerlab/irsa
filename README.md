## irsa
This reposetory contains two projects:
First, the IRSA project, which is a command line program to align theoretical IR/VCD to experimental spectra.
Second, the IRSA_GUI project, which provides a graphical user interface.


## Tutorial
The files for this tutorial can be found in the folder example. As a conformer sampler, we will use RDKit, and for the quantum mechanical frequency computations, we will use orca 4.2.0, since both programs are free of charge. Any other QM program (e.g. Gaussian) or conformer sampler (e.g. Omega) can also be used.

## First Considerations
For rigid and semi rigid compounds, the usage of the alignment algorithm should be straight-forward, and give satisfactory results. The results depend on the level of theory on which the quantum mechanical computations will be performed, and we can usually recommend BP86 functional with a triple zeta basis set (e.g. def2-tzvp, or cc-pVTZ) and a dispersion correction (e.g. Grimme's Dispersion correction with Becke Johnson damping). The most demanding step in the alignment procedure is the computation of the frequency spectra, since it requires the computation of the hessian matrix. For flexible compounds, this cost can quickly get quite demanding. To reduce the cost for flexible compounds, we thus recommend to lower the quality of the basis set (e.g. to def2-SVP). The first example we are going to discuss is the rigid compound Fenchone, others examples follow.

## Input files
Fenchone has two stereocenters, e.g., 2 Diastereomers could be distinguished using IR (the other isomers are entantiomers, i.e. they would require VCD to be detected). The smiles strings are:
1. [C@]12(C(=O)C([C@@H](C1)CC2)(C)C)C
2. [C@@]12(C(=O)C([C@@H](C1)CC2)(C)C)C

First, we generate the 3D input coordinates using the script 
```
import os
import subprocess
import sys
from rdkit import Chem
from rdkit.ML.Cluster import Butina
from rdkit.Chem import AllChem
from concurrent import futures
from rdkit import *
def get_conformer_energies(forcefield, mol,id):
    mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94')
    ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=id)
    energy = ff.CalcEnergy()
    return energy

fopen=open("smiles","r")
smiles=fopen.readline().split()[0]
mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)
ids = AllChem.EmbedMultipleConfs(mol,numConfs=5,pruneRmsThresh=0.001,randomSeed=42,numThreads=1,enforceChirality=True,useExpTorsionAnglePrefs=True,useBasicKnowledge=True) 
AllChem.MMFFOptimizeMoleculeConfs(mol,numThreads=1, mmffVariant='MMFF94')
mol_2 = Chem.RemoveHs(mol)
rmsmat = AllChem.GetConformerRMSMatrix(mol_2, prealigned=False)
num = mol.GetNumConformers()
rms_clusters = Butina.ClusterData(rmsmat, num, 1, isDistData=True, reordering=True)
for id2 in range(0,len(rms_clusters)):
    coord = mol.GetConformers()[id2].GetPositions()
    f = open(str(id2)+".xyz","w+")
    f.write(len(coord)+"\n"+str(id2)+"\n")
    for i in range(len(coord)):
        l = mol.GetAtomWidthIdx(i).GetSymbol()
        f.write(l+" "+str(coord[i][0])+" "+str(coord[i][1])+" "+str(coord[i][2])+"\n")
    f.close()conformer_rdkit.py
```















