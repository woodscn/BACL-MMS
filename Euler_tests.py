import os, sys

import numpy
import sympy
import thread

libpath = os.path.abspath('/Users/woodscn/')
sys.path.insert(0,libpath)
libpath = os.path.abspath('/home/woodscn/')
sys.path.insert(0,libpath)
from manufactured import Euler_UCS
xmin,xmax,ymin,ymax,zmin,zmax = 0,100,0,1,0,1
nx = 100
ny = 1
nz = 1
dxis = [1,1,1]
dt = 1./3.*dxis[0]
nt = 100
Euler_UCS = Euler_UCS.Euler_UCS(
    Euler_UCS.MASA_solution_full(
        ranges=[[xmin,xmax],[ymin,ymax],[zmin,zmax]],nxes=(nx,ny,nz),dxis=dxis,
        disc=True))
manufactured_source_function = Euler_UCS.balance_lambda_init()

t = sympy.Symbol('t')
xi = sympy.Symbol('xi')
eta = sympy.Symbol('eta')
zeta = sympy.Symbol('zeta')

def input_thread(L):
    raw_input()
    L.append(None)

f = open('euler_IMMS_terms.dat','w')
f.write(str(nx))
L = []
print "Press <enter> to end loop after current time step completes."
thread.start_new_thread(input_thread,(L,))
for indt in range(nt):
    for indx in range(nx):
        cell_ranges = [[t,indt*dt,(indt+1)*dt],
                       [xi,indx*dxis[0],(indx+1)*dxis[0]],
                       [eta,ymin,ymax],[zeta,zmin,zmax]]
        S_prime = Euler_UCS.balance_integrate(cell_ranges)
        f.write('\n'+', '.join([str(item) for item in (
            S_prime[0],S_prime[1],S_prime[2],S_prime[3],S_prime[4])]))
    f.write('\n')
    if L:
        break
f.write(str(indt))
f.close()

