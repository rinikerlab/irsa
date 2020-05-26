## irsa
This reposetory contains multiple scripts for the analysis of IR and VCD spectra:
First, the IRSA project, which is a command line program to align theoretical IR/VCD to experimental spectra. This is particular useful if one is interested in scanning automatically large data basis.
Second, the IRSA_GUI project, which provides a graphical user interface. This is particular useful if one is interested in an exact determination of the (absolute) stereochemistry.
Third, a toolbox, which purpose is to generate 
a) Conformers, b) perform baseline corrections for noisy spectra, c) generate input files for the IRSA/IRSA_GUI project.


## Tutorial
The files for this tutorial can be found in the folder example. As a conformer sampler, we will use RDKit, and for the quantum mechanical frequency computations, we will use orca 4.2.0, since both programs are free of charge. Any other QM program (e.g. Gaussian) or conformer sampler (e.g. Omega) can also be used.

## First Considerations
For rigid and semi rigid compounds, the usage of the alignment algorithm should be straight-forward, and should give satisfactory results. The results depend on the level of theory on which the quantum mechanical computations will be performed, and we can usually recommend BP86 functional with a triple zeta basis set (e.g. def2-tzvp, or cc-pVTZ) and a dispersion correction (e.g. Grimme's Dispersion correction with Becke Johnson damping, D3BJ). The most demanding step in the alignment procedure is the computation of the frequency spectra, since it requires the computation of the hessian matrix. For flexible compounds, this cost can quickly get quite demanding. To reduce the cost for flexible compounds, we thus recommend to lower the quality of the basis set (e.g. to def2-SVP). Obviously, this reduces the quality of the final spectrum obtained. The first example we are going to discuss is the rigid compound Fenchone, others examples follow.

## Generation of 3D Coordinates
Fenchone has two stereocenters, however, only one diastereomer makes chemically sence. The smile string is:
1. [C@]12(C(=O)C([C@@H](C1)CC2)(C)C)C

First, we generate the 3D input coordinates using the script conformer_rdkit.py, which creates possible 3D coordinates for each compound. It can be called via 
```
python conformer_rdkit.py "[C@]12(C(=O)C([C@@H](C1)CC2)(C)C)C"
```
and will create one conformer for the compound (since Fenchone is rigid).
We notice that the compounds are filtered via an heavy atom filtering, i.e. hydrogen atoms are not considered. For compounds, which have an hydroxyl group, we recommend to consider the Hydrogen atom.

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

smiles = sys.argv[1]
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
    f = open(str(id2)+"_unopt.xyz","w+")
    f.write(len(coord)+"\n"+str(id2)+"\n")
    for i in range(len(coord)):
        l = mol.GetAtomWidthIdx(i).GetSymbol()
        f.write(l+" "+str(coord[i][0])+" "+str(coord[i][1])+" "+str(coord[i][2])+"\n")
    f.close()
```
The script writes out .xyz files for each conformer found. We will use these coordinates as input files for the QM calculations.

## QM Computation

For the QM calculation, we will use the file ``0.inp``,
```
!RIJCOSX BP86 def2-SVP def2/J TightSCF TightOpt freq grid5 finalGrid6

* 0 1 xyzfile 0_unopt.xyz
```
, which uses the BP86/def2-SVP level of theory with the RIJCOSX approximation to perform the optimization and the frequency computation. Generally, dispersion correction (like D3BJ) is recommended. 
We will start the computation with
``` orca "0.inp > 0.out"```
, and the results will be stored in the 0.hess (Hessian), 0.engrad (gradient and energy) and 0.out (the output-log) files.
The computation will take approximate 2 hours using a single core.
We notice that most of the computations require multiple CPUs to be feasible. In such case you need to specifiy the number of processors avaible (via ```%pal nprocs N end```) and the memory available per processor (via ```%Mem N```) to speed up things considerable. 

## Generating the input files for the alignment algorithm and perform the alignment
The alignment algorithm reads pickle files stored in a specific format:
First, the (free) energy, saved as a numpy array of the format ```(number_of_conformers, 1)```.
Second, the frequency calculation output, saved as a numpy array of the format ```(number_of_conformers, number_of_peaks, 2)```. Here, the third axis saves the frequency of the normal-mode in [0], and the dipol or IR-intensity in [1].

The project comes with scripts, which convert the files to the desired format, for orca, it is ```convert_orca.py```.
```

```
The code goes through all ```.engrad``` and ```.hess``` files, and read the energies (quantum, in-vacuum), free energies and the frequencies and outputs the files energy.p, gibbs.p, freq.p


The alignment code can then be run by 
```python Main.py -e energy.p -f freq.p -o test```
and the code outputs the unaligned (testunshifted.out), aligned (testaligned.out), as well as the scores and the calculation options with the standard computation options. A list of computation options is provided at the end of this README.

## Using the GUI for alignment
For simple cases and screening purposes, the command line code is sufficient. However, for complex data, it might be desireable to pick peaks manually and analyse the spectrum more thorougly. To call the GUI, just call ```python Main_gui.py```














