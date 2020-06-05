## Version Beta 0.1

Please send requests and Bug reports to
Lennard.Boeselt@phys.chem.ethz.ch

If you use the algorithm or experimental data from this side and report your results, we would kindly ask you to cite

@article{Boeselt2019,
title={Determination of Absolute Stereochemistry of Flexible Molecules Using a Vibrational Circular Dichroism (VCD) Spectra Alignment Algorithm},
author={L. Boeselt and D. Sidler and T. Kittelmann and J. Stohner and D. Zindel and T. Wagner and S. Riniker},
journal={JCIM},
volume={59},
number={5},
pages={1826-1838},
year={2019},
doi = {10.1021/acs.jcim.8b00789}}  


@article{Boeselt2020,
title={Determination of Absolute Stereochemistry of Flexible Molecules Using a Vibrational Circular Dichroism (VCD) Spectra Alignment Algorithm},
author={L. Boeselt and R. Doetzer and S. Steiner and M. Stritzinger and S. Salzmann and S. Riniker},
journal={Determining the Regiochemistry and Relative Stereochemistry of Small and Druglike Molecules Using an Alignment Algorithm forvInfrared Spectra},
volume={xxx},
number={xxx},
pages={xxx},
year={2020},
doi = {10.1021/acs.analchem.0c01399}}  

and cite the paper where you got the data from.

If you use the provided Data or other scripts to analyze your IR spectra, we would kindly ask you to cite

@misc{rinikerlab2020,
title={Toolkit for IR analysis},
url={https://github.com/rinikerlab/irsa},
note={Software available from https://github.com/rinikerlab/irsa},
author={L. Boeselt and R. Doetzer and S. Steiner and M. Stritzinger and S. Salzmann and S. Riniker},
  year={2020},
}


## irsa
This reposetory contains multiple scripts for the analysis of IR and VCD spectra:
First, the IRSA project, which is a command line program to align theoretical IR/VCD to experimental spectra. This is particular useful if one is interested in scanning automatically large data basis.
Second, the IRSA_GUI project, which provides a graphical user interface. This is particular useful if one is interested in an exact determination of the (absolute) stereochemistry.
Third, a toolbox, which purpose is to 
a) generate Conformers, b) perform baseline corrections for noisy spectra, c) generate input files for the IRSA/IRSA_GUI project.


## Tutorial
The files for this tutorial can be found in the folder example. As a conformer sampler, we will use RDKit, and for the quantum mechanical frequency computations, we will use orca 4.2.0, since both programs are free of charge. Any other QM program (e.g. Gaussian) or conformer sampler (e.g. Omega) can also be used. Python scripts to sample conformers with RDKit and with Omega are available in the Toolbox. Scripts for converting orca and gaussian output files to the necessary pickle files are also available.

## First Considerations
For rigid and semi rigid compounds, the usage of the alignment algorithm should be straight-forward. From the experimental side, we can recommend to perform GC-IR experiments to remove the influence from the solvent, which often interacts with the solute. If this is not possible, or if the GC-IR spectrum is very noisy, IR/VCD experiments can be performed in solvent. Here, it is often desirable to try different solvents (Chloroform, DMSO). Further, we have the experience that especially the region below 1000 wavenumbers strongly differs between isomers. Therefore, it is desirable to aim for a high resolution in this region. **Wavenumbers above 1500 are usually not suitable for determining the stereochemistry**.  From the computational side, the results depend on the level of theory on which the quantum mechanical computations will be performed, and we can usually recommend the BP86 functional with the cc-pVTZ basis set for rigid conformers. This combination is known to give a favourable error cancellation. However, the basis set is quite expensive, and if the computational cost is a limiting factor, changing the basis set to def2-tzvp or def2-SVP might be considered. For flexible compounds, we also need to achieve chemical accuracy, therefore switching to a hybrid functional might become necessary. The most demanding step in the alignment procedure is the computation of the frequency spectra, since it requires the computation of the hessian matrix. For flexible compounds, this cost can quickly get quite demanding. Obviously, this reduces the quality of the final spectrum obtained. The first example we are going to discuss is the rigid compound Fenchone, others examples follow.

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
!RI BP86 def2-tzvp def2/J TightSCF TightOpt D3BJ freq grid5 finalGrid6
%pal nprocs 12
end

* 0 1 xyzfile 0_unopt.xyz
```
, which uses the BP86/def2-tzvp level of theory with the RIJCOSX approximation to perform the optimization and the frequency computation. Generally, dispersion correction (like D3BJ) is recommended. 
We will start the computation with
``` bsub -n 12 "orca 0.inp > 0.out"```
, and the results will be stored in the 0.hess (Hessian), 0.engrad (gradient and energy) and 0.out (the output-log) files.
However, the command line might vary depending on your cluster architecture. If you only have a single core available, I would recommend to reduce the level of theory, e.g.,

```
!RI BP86 def2-SVP def2/J TightSCF TightOpt D3BJ freq grid5 finalGrid6
%pal nprocs 12
end

* 0 1 xyzfile 0_unopt.xyz
```

We notice that most of the computations require multiple CPUs to be feasible. In such case you need to specifiy the number of processors avaible (via ```%pal nprocs N end```) and the memory available per processor (via ```%Mem N```) to speed up things considerable. 

## Generating the input files for the alignment algorithm and perform the alignment
The alignment algorithm reads pickle files stored in a specific format:
First, the (free) energy, saved as a numpy array of the format ```(number_of_conformers, 1)```.
Second, the frequency calculation output, saved as a numpy array of the format ```(number_of_conformers, number_of_peaks, 2)```. Here, the third axis saves the frequency of the normal-mode in [0], and the dipol or IR-intensity in [1].

The project comes with scripts, which convert the files to the desired format, for orca (and gaussian), it is ```convert_orca.py``` and ```convert_gaussian.py```.
The python scripts take two inputs: The number of atoms and the directory where to find the output files, e.g.,
```
python convert_orca.py -n 27 -d /DIR/TO/YOUR/CALC/
```
The code goes through all ```.out``` files, and read the energies (quantum, in-vacuum), free energies and the frequencies and outputs the files energy.p, gibbs.p, freq.p, which are necessary for the alignment


## Performing the alignment with the command line version

The alignment code can then be run by 
```python Main.py -e energy.p -f freq.p -exp IR.txt -lb 1000 -hb 1500 -o test -mu 0.99 -s1 0.1-s2 0.01```
and the code outputs the unaligned (testunshifted.out), aligned (testaligned.out), as well as the scores (test.txt) and the calculation options with the standard computation options. A list of computation options is provided at the end of this README.



## Performing the alignment with the Graphical User Interface


To use the GUI instead, one needs to call the program
```python Main_GUI.py```.
[...] OUR CLUSTER IS DOWN, THUS THIS TUTORIAL NEEDS APPROX ONE MORE WEEK TO BE FINISHED [...]







## COMPUTATION OPTIONS
