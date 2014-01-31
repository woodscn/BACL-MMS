import sympy
from functools import partial

from heat_equation import MASA_solution, HeatEquation, MASA_source_lambda
from Euler_UCS import Euler_UCS, unsteady_Euler

class Error(Exception):
    pass


def RD_test(t,x,y,z):
    return t**2*x*y*z

if __name__=="__main__":
#    print recursive_derivative(RD_test,(1,0,0,1))-2.0
    t = sympy.Symbol('t')
    x = sympy.Symbol('x')
    y = sympy.Symbol('y')
    z = sympy.Symbol('z')
    xi = sympy.Symbol('xi')
    eta = sympy.Symbol('eta')
    zeta = sympy.Symbol('zeta')
    args = (.01,.01,.01,.01)
    dxes = [.001 for arg in args]
#    kwargs = {
#        'Ax':1,'At':1,'By':.5,'Bt':-.25,'Cz':.7,'Ct':0,'Dt':.1,'rho':1,'cp':1,'k':1}
#    eqn = HeatEquation(MASA_solution(**kwargs))
#    print abs(recursive_derivative(
#            lambda x0,x1,x2,x3:eqn.balance_integrate(
#                ((t,0,x0),(x,0,x1),(y,0,x2),(z,0,x3))),args,dxes,order=5) - 
#              MASA_source_lambda(**kwargs)(*args))
    eqn2 = Euler_UCS(unsteady_Euler('two_shock'))
    tmin, tmax = 0.0001, .1
    xmin, xmax = -1, 1
    ymin, ymax = -1, 1
    zmin, zmax = -1, 1

# # Debugging stuff here.
#     import numpy
#     import pylab 
#     import functools
#     from sympy.utilities.lambdify import lambdify
#     x=pylab.linspace(xmin,xmax)
#     y=pylab.linspace(ymin,ymax)
#     z=pylab.linspace(zmin,zmax)
#     tarray=pylab.linspace(tmin,tmax)
#     coord1, coord2, coord3, coord4 = (xi,x), (eta,y), (t,tarray), (zeta,z)
#     grid1, grid2 = numpy.meshgrid(coord1[1], coord2[1])
#     evaluate_sol = lambdify((coord3[0],coord4[0],coord1[0],coord2[0]),
#                             eqn2.sol[1])
#     partial_sol = functools.partial(evaluate_sol,.05,1)
#     grid_sol = numpy.vectorize(partial_sol)(grid1,grid2)
#     junk = pylab.contourf(grid1,grid2,grid_sol)
#     pylab.colorbar(junk)
#     pylab.show()
#     import pdb;pdb.set_trace()

    print eqn2.balance_integrate(
                ((t,tmin,tmax)
                 ,(xi,xmin,xmax)
                 ,(eta,ymin,ymax)
#                 ,(zeta,zmin,zmax)
                 ))
    print "done"
              

