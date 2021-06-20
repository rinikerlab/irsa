## Version Beta 0.2

We implemented the algorithm, such that it works with the convoluted theoretical spectra, instead of the output of the QM software packages. It is recommended to use this code instead of the previous implementation. It follows a short tutorial.

## Test case
![alt text](https://github.com/rinikerlab/irsa/blob/master/diastereomeric_pair.png)
Consider the diastereomeric pair 1a and 2a, where we measured the experimental IR and VCD spectrum to determine the relative stereochemistry and the absolute stereochemistry. We will refer to the enantiomer of 1a as 1b, and to the enantiomer of 2a as 2b.
First, we generate the conformational ensemble using RDKit.
The isomeric SMILES string for compound 1a is
O[C@]1([H])[C@](C)(C2(C)C)CC[C@H]2C1
and the isomeric SMILES string for compound 2a is
[H][C@]1(O)[C@](C2(C)C)(CC[C@H]2C1)C
. Save the string of compound 1a to a file called ```SMILES```.
The code snipped
```
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
    try:
        rdkit.Chem.rdmolfiles.MolToXYZFile(solute_mol, str(i)+".xyz", confId=i)
    except:
        print("Exhausted")
        break      
```
will generate a conformational ensemble of 1a and 2a. Simply copy & paste the code to a file (we will name it ```create.py```), and call 
```
python create.py
```
. This will generate n xyz-files, one for each conformer found.
For each conformer, prepare a computational input file. We will use gaussian, version 09D.
For compound 1a, this may look like
```
%mem=16000MB
%NProcShared=16
#B3LYP/Gen scf=tight freq=VCD Opt=Tight Integral(Grid=SuperFine) EmpiricalDispersion=GD3

0

0 1
O      1.902578   -0.390373   -1.506707
C      1.143131   -1.021854   -0.523900
C     -0.299644   -0.559541   -0.576112
C     -0.922356   -0.577159   -1.894603
C     -0.187104    0.797228    0.083466
C      0.707782    1.743241   -0.639092
C     -1.504518    1.429500    0.363166
C     -1.036367   -1.340693    0.502571
C     -0.563181   -0.725133    1.825531
C      0.434573    0.290022    1.393955
C      1.628243   -0.467675    0.822702
H      2.828656   -0.306770   -1.155892
H      1.196109   -2.113809   -0.564316
H     -1.945529   -1.048327   -1.818866
H     -0.353505   -1.291899   -2.559489
H     -1.021334    0.372263   -2.425153
H      0.322476    2.797734   -0.513075
H      1.752207    1.720400   -0.326854
H      0.603401    1.541416   -1.733088
H     -2.051490    1.667494   -0.555298
H     -1.323200    2.382327    0.892540
H     -2.176397    0.833684    1.000923
H     -0.660822   -2.397316    0.406226
H     -2.119211   -1.313686    0.393149
H     -0.006732   -1.527789    2.380344
H     -1.380448   -0.367860    2.442272
H      0.685277    1.067176    2.083587
H      1.828339   -1.340072    1.479451
H      2.519066    0.147472    0.722558

@basis.gbs

```
where basis.gbs refers to the def2-TZVP basis set and can be extracted from the folder tutorial. The computation will result in .log files, which contain all information from the computation.
Next, we extract the IR and the VCD spectrum from the calculation. Simply place the script ```IR_and_VCD_g09.py``` in the folder of your log files, type ```python IR_and_VCD_g09.py```, which will result in the ```IR_theo.txt``` and ```vcd_theo.txt```. Place the experimental IR and VCD spectra together with the theoretical spectra in one folder, and call ```python pick_peaks_exp_IR_VCD.py```. This will pick the peaks of the experimental and theoretical IR and VCD spectra by a simple distance criterium, and will create the files ```peaks_mod_theo.txt```,```peaks_mod_exp.txt```, ```peaks_mod_theo_vcd.txt``` and ```peaks_mod_exp_vcd.txt```. In each file, the first column refers to the normalized intensity, the second column to the x-axis, and the third column refers to, whether the data point belongs to the IR spectrum (0) or to the VCD spectrum (1). You can create these files also on your own, when you have noisy data and the peak picking algorithm does not perform appropriately.
Next, we are going to align the spectra to each other. For this, consider the file
```Settings.py```, which reads
```
class Settings:
    def __init__(self):
        self.directory='1a'
        self.use_vcd=True

        self.shift_factor=0.98
        self.cutoff_absolute=False
        self.cutoff=0.015
        self.sigma_1=0.192
        self.sigma_2=0.1
        self.x_min=1000
        self.x_max=1500
    def get(self):
        return self.shift_factor, self.cutoff_absolute, self.cutoff, self.sigma_1, self.sigma_2, self.use_vcd, self.x_min, self.x_max
    def get_directory(self):
        return self.directory
```
.
Most importantly is the variable ```self.directory```, which refers to the total path of the folder, where the files ```IR_theo.txt```,```vcd_theo.txt```, ```peaks_mod_theo.txt```,```peaks_mod_exp.txt```, ```peaks_mod_theo_vcd.txt```, ```peaks_mod_exp_vcd.txt``` and the experimental spectra ```VCD.txt``` and ```IR.txt``` are saved. Next, call ```python Main.py```.
This will output the png files ```1a_1_0.98_0.015_vcd__1.png```,```1a_1_0.98_0.015_vcd__-1.png```, and ```1b_0_0.98_0.015_vcd_.png```. The first spectrum belong to the aligned vcd spectrum of the provided diastereomer, the second file to the aligned vcd spectrum of the enatiomer, and the third file to the aligned IR spectrum. Black refers always to the experimental spectrum, red to the aligned spectrum, and orange to the unshifted, but by the shift-factor scaled spectrum. Further, it will output information of the alignment as
```
VCD SPECTRUM AVAILABLE
CALL ALIGNMENT
Initialization ND
Initialization SUCCESSFUL
START ALIGNMENT
ALIGNMENT SUCCEESFUL
returnvalue -0.44260629081427527
READ IN EXPERIMENT
READ IN THEORETICAL SPECTRUM
Fontconfig warning: ignoring UTF-8: not a valid region tag
-1 Alignment Score 0.44260629081427527
-1 IR pearson number 0.2479916761927686
-1 IR spearman number 0.5056763587054348
-1 VCD pearson number 0.088569617048474
-1 VCD spearman number 0.1679884799539198
VCD SPECTRUM AVAILABLE
CALL ALIGNMENT
Initialization ND
Initialization SUCCESSFUL
START ALIGNMENT
ALIGNMENT SUCCEESFUL
returnvalue -0.4421465263852314
READ IN EXPERIMENT
READ IN THEORETICAL SPECTRUM
1 Alignment Score 0.4421465263852314
1 IR pearson number 0.20554358106704632
1 IR spearman number 0.37500399601598405
1 VCD pearson number 0.13037652458387428
1 VCD spearman number 0.12164679058716235
sucess
```
, where the pearson number and spearman number refers to overlap metrics between the experimental and theoretical IR spectrum, as well as the experimental and theoretical VCD spectrum. The returnvalue refers to the score of the alignment. The -1 refers to the enantiomer of the provided compound (in this case, 1b).
The aligned spectra will look like
![alt text](https://github.com/rinikerlab/irsa/blob/master/1a_IR.png)
![alt text](https://github.com/rinikerlab/irsa/blob/master/1a_VCD.png)
![alt text](https://github.com/rinikerlab/irsa/blob/master/1b_VCD.png).
Repeating this procedure with compound 2a will result in 
```
VCD SPECTRUM AVAILABLE
CALL ALIGNMENT
Initialization ND
Initialization SUCCESSFUL
START ALIGNMENT
ALIGNMENT SUCCEESFUL
returnvalue -0.7702852131065266
READ IN EXPERIMENT
READ IN THEORETICAL SPECTRUM
Fontconfig warning: ignoring UTF-8: not a valid region tag
-1 Alignment Score 0.7702852131065266
-1 IR pearson number 0.8167838895785107
-1 IR spearman number 0.594695114780459
-1 VCD pearson number 0.820094471868076
-1 VCD spearman number 0.6286756027024107
VCD SPECTRUM AVAILABLE
CALL ALIGNMENT
Initialization ND
Initialization SUCCESSFUL
START ALIGNMENT
ALIGNMENT SUCCEESFUL
returnvalue -0.743434233344465
READ IN EXPERIMENT
READ IN THEORETICAL SPECTRUM
1 Alignment Score 0.743434233344465
1 IR pearson number 0.815046461294488
1 IR spearman number 0.5928068832275329
1 VCD pearson number -0.7752980148118074
1 VCD spearman number -0.5317340309361237
sucess
```
and the spectra
![alt text](https://github.com/rinikerlab/irsa/blob/master/1a_IR.png)
![alt text](https://github.com/rinikerlab/irsa/blob/master/1a_VCD.png)
![alt text](https://github.com/rinikerlab/irsa/blob/master/1b_VCD.png).

From the overlap metrics and the figures, it is clear that compound 2b is the compound which was measured.






 


## Data for the submission F. Pultar, JACS, 2021, submitted, can be found in the directory "example"

## Version Beta 0.1
We fixed a part in the Boltzmann weighting
This side is still under construction
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

Next, we perform the alignment of (+)-Borneol via the graphical user interface.
For this, you first need to install Qt5 for python. Installation instructions can be found on the offical website.

Go to the folder GUI and start the code via ```python Gui.py```.

A window will open, which will look like the pictures at the end in this chapter.
This window is divided into 3 subpanels (experimental spectra, theoretical spectra and alignment).
First, read in the experimental IR spectrum (DATA/borneol_IR.txt) and the experimental VCD spectraum (DATA/borneol_VCD.txt)
Next, click the button normalize to normalize the spectra to intensity 1.
Then, click the Button automatic peak selection to obtain the peaks.

Next, you will need in to load the theoretical data. For this, click on LOAD E (example_borneol/energy.p) then LOAD IR (example/freq.p) and LOAD VCD (example_borneol/vcd.p)
You can set the Lorentzband width to 12 to have an spectrum which agrees better with experiment.

Next, press change the value mu to 1.0135 (the spectrum was computed with BP86/cc-pVTZ, where the shift factor of 1.0135 is an appropiated value). Next, press align, and you results should look like the one in the pictures provided below.

The textfield to the right provides the numeric results obtained.
The score (more negative is better!), the pearson coefficient of the aligned IR and the pearson coefficient of the aligned VCD spectrum. If you only want to align IR spectra, you can just leave out the VCD spectrum.




![alt text](https://github.com/rinikerlab/irsa/blob/master/gui1.png)


![alt text](https://github.com/rinikerlab/irsa/blob/master/gui2.png)



## Computational options
