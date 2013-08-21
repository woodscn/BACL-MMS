import sympy
from functools import partial

from heat_equation import MASA_solution, HeatEquation, MASA_source_lambda
from Euler_UCS import Euler_UCS, MASA_solution_E, unsteady_Euler

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
    kwargs = {
        'Ax':1,'At':1,'By':.5,'Bt':-.25,'Cz':.7,'Ct':0,'Dt':.1,'rho':1,'cp':1,'k':1}
    eqn = HeatEquation(MASA_solution(**kwargs))
    print abs(recursive_derivative(
            lambda x0,x1,x2,x3:eqn.balance_integrate(
                ((t,0,x0),(x,0,x1),(y,0,x2),(z,0,x3))),args,dxes,order=5) - 
              MASA_source_lambda(**kwargs)(*args))
    import pdb;pdb.set_trace()
    eqn2 = Euler_UCS(unsteady_Euler('normal'))
    print "got this far"
    print eqn2.balance_integrate(
                ((t,0.9,1.),(xi,-1.,1.),(eta,-1.,1.)))#,(zeta,-1.,1.)))
    print "done"
              

