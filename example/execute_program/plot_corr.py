import numpy as np
import matplotlib.pyplot as py

py.rc('font',family='serif')
py.rc('xtick',labelsize='x-small')
py.rc('ytick',labelsize='x-small')
#1.0135
fig = py.figure(figsize=(4,3.0))
ax = fig.add_subplot(1,1,1)

#mu = np.arange(0.990,1.040,0.005)[:-1]
mu = np.arange(1.000,1.014,0.001)#[:-1]
print(mu)

pearson_I = np.loadtxt("I_unshifted.txt",usecols=(1,))[::2]
np.savetxt("I_unshifted_pearson_500_1000.txt",np.concatenate((mu[:,np.newaxis],pearson_I[:,np.newaxis]),axis=1))

pearson_II = np.loadtxt("II_unshifted.txt",usecols=(1,))[::2]
np.savetxt("II_unshifted_pearson_500_1000.txt",np.concatenate((mu[:,np.newaxis],pearson_II[:,np.newaxis]),axis=1))

spearman_I = np.loadtxt("I_unshifted.txt",usecols=(1,))[1::2]
np.savetxt("I_unshifted_spearman_500_1000.txt",np.concatenate((mu[:,np.newaxis],spearman_I[:,np.newaxis]),axis=1))

spearman_II = np.loadtxt("II_unshifted.txt",usecols=(1,))[1::2]
np.savetxt("II_unshifted_spearman_500_1000.txt",np.concatenate((mu[:,np.newaxis],spearman_II[:,np.newaxis]),axis=1))

ax.plot(mu,pearson_I,color='blue',label="I,Pearson")
ax.plot(mu,pearson_II,color='green',label="II,Pearson")
ax.plot(mu,spearman_I,'--',color='blue',label="I,Spearman")
ax.plot(mu,spearman_II,'--',color='green',label="II,Spearman")
ax.set_ylim(-0.1,1.0)
ax.set_xlabel('$\mu$',labelpad=-2)
ax.set_ylabel('corr / a.u.')
py.legend(loc=2,prop={'size':6})
fig.savefig('scores_500_1000.png')
py.show()
#
fig = py.figure(figsize=(4,3.0))
ax = fig.add_subplot(1,1,1)
mu = np.arange(1.000,1.014,0.001)#[:-1]

pearson_IRSA_I = np.loadtxt("I_shifted.txt",usecols=(1,))[::2]
np.savetxt("I_IRSA_pearson_500_1000.txt",np.concatenate((mu[:,np.newaxis],pearson_IRSA_I[:,np.newaxis]),axis=1))

pearson_IRSA_II = np.loadtxt("II_shifted.txt",usecols=(1,))[::2]
np.savetxt("II_IRSA_pearson_500_1000.txt",np.concatenate((mu[:,np.newaxis],pearson_IRSA_II[:,np.newaxis]),axis=1))

spearman_IRSA_I = np.loadtxt("I_shifted.txt",usecols=(1,))[1::2]
np.savetxt("I_IRSA_spearman_500_1000.txt",np.concatenate((mu[:,np.newaxis],spearman_IRSA_I[:,np.newaxis]),axis=1))

spearman_IRSA_II = np.loadtxt("II_shifted.txt",usecols=(1,))[1::2]
np.savetxt("II_IRSA_spearman_500_1000.txt",np.concatenate((mu[:,np.newaxis],spearman_IRSA_II[:,np.newaxis]),axis=1))

ax.plot(mu,pearson_IRSA_I,color='blue',label="I,IRSA-Pearson")
ax.plot(mu,pearson_IRSA_II,color='green',label="II,IRSA-Pearson")
ax.plot(mu,spearman_IRSA_I,"--",color='blue',label="I,IRSA-Spearman")
ax.plot(mu,spearman_IRSA_II,"--",color='green',label="II,IRSA-Spearman")
ax.set_ylim(-0.1,1.0)
ax.set_xlabel('$\mu$',labelpad=-2)
ax.set_ylabel('corr / a.u.')
py.legend(loc=2,prop={'size':6})
fig.savefig('scores_IRSA_500_1000.png')
py.show()
fig = py.figure(figsize=(4,3.0))
ax = fig.add_subplot(1,1,1)
mu = np.arange(1.000,1.014,0.001)#[:-1]
pearson_I = np.loadtxt("I_500_1800.txt",usecols=(1,))[::2]
np.savetxt("I_unshifted_pearson_500_1800.txt",np.concatenate((mu[:,np.newaxis],pearson_I[:,np.newaxis]),axis=1))

pearson_II = np.loadtxt("II_500_1800.txt",usecols=(1,))[::2]
np.savetxt("II_unshifted_pearson_500_1800.txt",np.concatenate((mu[:,np.newaxis],pearson_II[:,np.newaxis]),axis=1))

spearman_I = np.loadtxt("I_500_1800.txt",usecols=(1,))[1::2]
np.savetxt("I_unshifted_spearman_500_1800.txt",np.concatenate((mu[:,np.newaxis],spearman_I[:,np.newaxis]),axis=1))

spearman_II = np.loadtxt("II_500_1800.txt",usecols=(1,))[1::2]
np.savetxt("II_unshifted_spearman_500_1800.txt",np.concatenate((mu[:,np.newaxis],spearman_II[:,np.newaxis]),axis=1))


ax.plot(mu,pearson_I,color='blue',label="I,Pearson")
ax.plot(mu,pearson_II,color='green',label="II,Pearson")
ax.plot(mu,spearman_I,'--',color='blue',label="I,Spearman")
ax.plot(mu,spearman_II,'--',color='green',label="II,Spearman")
ax.set_ylim(-0.1,1.0)
ax.set_xlabel('$\mu$',labelpad=-2)
ax.set_ylabel('corr / a.u.')
py.legend(loc=2,prop={'size':6})
fig.savefig('scores_500_1800.png')
py.show()
fig = py.figure(figsize=(4,3.0))
ax = fig.add_subplot(1,1,1)
mu = np.arange(1.000,1.014,0.001)#[:-1]
##500 1800 IRSA
pearson_IRSA_I = np.loadtxt("I_500_1800_shifted.txt",usecols=(1,))[::2]#[::3]
np.savetxt("I_IRSA_pearson_500_1800.txt",np.concatenate((mu[:,np.newaxis],pearson_IRSA_I[:,np.newaxis]),axis=1))


pearson_IRSA_II = np.loadtxt("II_500_1800_shifted.txt",usecols=(1,))[::2]
np.savetxt("II_IRSA_pearson_500_1800.txt",np.concatenate((mu[:,np.newaxis],pearson_IRSA_II[:,np.newaxis]),axis=1))

spearman_IRSA_I = np.loadtxt("I_500_1800_shifted.txt",usecols=(1,))[1::2]#[1::3]
np.savetxt("I_IRSA_spearman_500_1800.txt",np.concatenate((mu[:,np.newaxis],spearman_IRSA_I[:,np.newaxis]),axis=1))

spearman_IRSA_II = np.loadtxt("II_500_1800_shifted.txt",usecols=(1,))[1::2]
np.savetxt("II_IRSA_spearman_500_1800.txt",np.concatenate((mu[:,np.newaxis],spearman_IRSA_II[:,np.newaxis]),axis=1))
ax.plot(mu,pearson_IRSA_I,color='blue',label="I,IRSA-Pearson")
ax.plot(mu,pearson_IRSA_II,color='green',label="II,IRSA-Pearson")
ax.plot(mu,spearman_IRSA_I,"--",color='blue',label="I,IRSA-Spearman")
ax.plot(mu,spearman_IRSA_II,"--",color='green',label="II,IRSA-Spearman")
ax.set_ylim(-0.1,1.0)
ax.set_xlabel('$\mu$',labelpad=-2)
ax.set_ylabel('corr / a.u.')
py.legend(loc=2,prop={'size':6})
fig.savefig('scores_IRSA_500_1800.png')
py.show()


