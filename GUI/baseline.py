import numpy as np
import scipy
import matplotlib.pyplot as py
a = np.loadtxt("IR_absorbance.txt",usecols=(0,1,))
b_1 = 500
b_2 = 900
vals =  (a[:,0] >= b_1) & (b_2 >=a[:,0])

p = np.polyfit(a[vals,0],a[vals,1]-np.min(a[vals,1]),0)
p = np.poly1d(p)
arr = p(a[vals,0])
a[vals,1] -= arr
a[vals,1] -= np.min(a[vals,1])
py.plot(a[vals,0],a[vals,1])
py.plot(a[vals,0],arr)
py.show()
np.savetxt("IR_absorbance_baseline_"+str(b_1)+"_"+str(b_2)+".txt",a[vals])
