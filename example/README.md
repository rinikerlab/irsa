Data for the submitted Article <br/>
F. Pultar, JACS, 2021, submitted.
<br/>
The structure is as follows:<br/>
Directory OMEGA:<br/>
  Sampling of conformers<br/>
DFT:<br/>
  Input files and output files from the DFT calculations<br/>
SCRIPTS:<br/>
  python script, to create the necessary pickle files<br/>
<br/>
<br/>
Notice that we did not know the correct structure at this stage. Hence, the archives are named 0 (isomer epi-I in the publication) and 1 (isomer I, correct isomer).

The essential part is in the directory execute_program:<br>
Here, the baseline-corrected experimental spectrum is provided, as well as the hand picked peaks. <br>
The script run_constant.sh reads the provided pickle files, and runs the computation and the alignment with the provided parameters (see run_constant file).
Notice that we focus on the found minimum structure, because the conformers are well separated in energy. You can pickle the files yourself using the script convert_freq_ZPE.py provided. This will alter the results a little bit, since the complete ensemble will be used. <br>
