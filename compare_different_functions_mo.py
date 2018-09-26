import numpy as np
import matplotlib.pyplot as plt

def file_len(fname):
    with open(fname) as f:
       for i, l in enumerate(f):
           pass
    return i + 1

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=22)
params = {'legend.fontsize': 20}
plt.rcParams.update(params)

frho1 = open('test-rhoc-mo', 'r')
nrho1 = file_len('test-rhoc-mo')-2 # minus 2, because first lines are comments
frho2 = open('test-rhoc-mo-parsec-cfac2', 'r')
nrho2 = file_len('test-rhoc-mo-parsec-cfac2')-2

rho1 = np.zeros((nrho1,2))
rho2 = np.zeros((nrho2,2))
#skip first two lines
for i in range(2):
    line = frho1.readline()
    line = frho2.readline()

for i in range(nrho1):
    line = (frho1.readline()).split()
    rho1[i,0] = float(line[1])
    rho1[i,1] = float(line[2])

for i in range(nrho2):
    line = (frho2.readline()).split()
    rho2[i,0] = float(line[1])
    rho2[i,1] = float(line[2])

t = np.arange(0,1200,10)
s = np.zeros(len(t))
plt.xlabel(r'$|G|^2/2$ (Rydberg)')
plt.ylabel(r'$|\rho^{partial}_{core}(G)|$ (a.u.)$^{-3}$')
plt.plot(rho1[:,0],abs(rho1[:,1])/4.2875,'r.',markersize=2.5)
plt.plot(rho2[:,0],abs(rho2[:,1])/4.2875,'b.',markersize=2.5)
plt.plot(-1,100.0,'r.',markersize=6,label=r'$\rho^{partial}_{core}(r<r_{core}) = \frac{A\sin(Br)}{r}$')
plt.plot(-1,100.0,'b.',markersize=6,label=r'$\rho^{partial}_{core}(r<r_{core}) = \exp(c_0+c_1r^2+c_2r^4)$')
plt.legend(loc='upper right', fontsize='small')
plt.plot(t,s,'k--')

plt.xlim(400,1000)
plt.ylim(-0.1e-8,0.5e-7)
plt.tight_layout()
#plt.show()
#plt.savefig("compare_different_functions_mo.pdf")
plt.savefig("compare_different_functions_mo-cfac2.png",dpi=300)
#plt.savefig("compare_different_functions_mo.eps")
