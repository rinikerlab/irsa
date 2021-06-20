%mem=16000MB
%NProcShared=24
#B3LYP/Gen scf=tight freq=VCD Opt=Tight Integral(Grid=SuperFine) EmpiricalDispersion=GD3

0

0 1
@/cluster/work/igc/boeseltl/NEW_ISOMERS/Borneol_2/0.xyz

@/cluster/work/igc/boeseltl/NEW_ISOMERS/Borneol/basis.gbs
    
    
