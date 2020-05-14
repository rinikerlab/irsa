## irsa
This reposetory contains two projects:
First, the IRSA project, which is a command line program to align theoretical IR/VCD to experimental spectra.
Second, the IRSA_GUI project, which provides a graphical user interface.


## Tutorial
The files for this tutorial can be found in the folder example. As a conformer sampler, we will use RDKit, and for the quantum mechanical frequency computations, we will use orca 4.2.0, since both programs are free of charge. Any other QM program (e.g. Gaussian) or conformer sampler (e.g. Omega) can also be used.

## First Considerations
For rigid and semi rigid compounds, the usage of the alignment algorithm should be straight-forward, and give satisfactory results. The results depend on the level of theory on which the quantum mechanical computations will be performed, and we can usually recommend BP86 functional with a triple zeta basis set (e.g. def2-tzvp, or cc-pVTZ) and a dispersion correction (e.g. Grimme's Dispersion correction with Becke Johnson damping). The most demanding step in the alignment procedure is the computation of the frequency spectra, since it requires the computation of the hessian matrix. For flexible compounds, this cost can quickly get quite demanding. To reduce the cost for flexible compounds, we thus recommend to lower the quality of the basis set (e.g. to def2-SVP). The first example we are going to discuss is the rigid compound Fenchone, others examples follow

## Conformer generation










