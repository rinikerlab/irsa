import sys
import random
import numpy as np
import argparse
import Reader
import Spectrum
import Algorithm
values = {"w":6,"lb":1000,"hb":1600,"mu":0.98,"s0":0.1,"s1":0.002,"c":999,"T":300}
load = {"load_exp_ir":"","load_exp_vcd":"","load_theo_ir":"","load_theo_exp":""}
exp_ir = 0
set_IR = False
set_VCD = False
T = 300
Dipol = 0
