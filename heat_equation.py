from itertools import product

import sympy
from sympy.utilities.lambdify import lambdify
from sympy import sin, cos

import numpy
import numpy.random

from base_equation import SympyEquation
from recursive_derivative import recursive_derivative

Zero = sympy.singleton.S.Zero
t = sympy.Symbol('t')
x = sympy.Symbol('x')
y = sympy.Symbol('y')
z = sympy.Symbol('z')
xi= sympy.Symbol('xi')

class HeatEquation(SympyEquation):
    def setup(self,**kwargs):
        self.rho = kwargs['rho']
        self.cp = kwargs['cp']
        self.k = kwargs['k']
        self.fluxes = [sympy.Matrix([self.rho*self.cp*self.sol[0]]),
                       sympy.Matrix(
                [-self.k*sympy.diff(self.sol[0],self.vars_[1])]),
                       sympy.Matrix(
                [-self.k*sympy.diff(self.sol[0],self.vars_[2])]),
                       sympy.Matrix(
                [-self.k*sympy.diff(self.sol[0],self.vars_[3])])]
        self.source = sympy.Matrix([Zero])

def MASA_solution(Ax,At,By,Bt,Cz,Ct,Dt,rho,cp,k):
    return {'vars':[t,x,y,z],'eqn_kwargs':{'rho':rho,'cp':cp,'k':k},
                'sol':[sympy.cos(Ax*x+At*t)*sympy.cos(By*y+Bt*t)*
                       sympy.cos(Cz*z+Ct*t)*sympy.cos(Dt*t)],
                'discontinuities':[]}

def MASA_source(Ax,At,By,Bt,Cz,Ct,Dt,rho,cp,k):
    out = (
        -(sin(Ax*x+At*t)*cos(By*y+Bt*t)*cos(Cz*z+Ct*t)*cos(Dt*t)*At+
          cos(Ax*x+At*t)*sin(By*y+Bt*t)*cos(Cz*z+Ct*t)*cos(Dt*t)*Bt+
          cos(Ax*x+At*t)*cos(By*y+Bt*t)*sin(Cz*z+Ct*t)*cos(Dt*t)*Ct+
          cos(Ax*x+At*t)*cos(By*y+Bt*t)*cos(Cz*z+Ct*t)*sin(Dt*t)*Dt)*rho*cp+
         (Ax**2+By**2+Cz**2)*cos(Ax*x+At*t)*cos(By*y+Bt*t)*cos(Cz*z+Ct*t)*
         cos(Dt*t)*k
         )
    return out

def MASA_source_lambda(**kwargs):
    return lambdify((t,x,y,z),MASA_source(**kwargs))

def heat_exact_sol(space,n):
    # These solutions taken from:
    # http://eqworld.ipmnet.ru/en/solutions/lpde/lpde101.pdf
    A, B, C, mu = 1,1,1,1
    a = 1
    out= [A*space+B,
          A*(space**2+2*a*t)+B,
          A*(space**3+6*a*t*space)+B,
          A*(space**4+12*a*t*space**2+12*a**2*t**2)+B,
          A*sympy.exp(a*mu**2*t+mu*space)+B,
          A*sympy.exp(a*mu**2*t-mu*space)+B,
          A*sympy.exp(-a*mu**2*t)*sympy.cos(mu*space+B)+C,
          A*sympy.exp(-mu*space)*sympy.cos(mu*space-2*a*mu**2*t+B)+C]
    return out[n]

def heat_exact_tests(ranges,nx):
    angle_list = [range_*(.5*sympy.pi/(nx-1)) for range_ in range(nx)]
    spaces = [sympy.cos(angle)*sympy.cos(angle)*x
              -sympy.cos(angle)*sympy.sin(angle)*y
              +sympy.sin(angle)*z for angle in angle_list]
    return [[heat_exact_sol(xi,n),angle_list,heat_exact_array(spaces,n,ranges)]
            for n in range(8)]
 

def heat_exact_array(spaces,n,ranges):
    return [heat_residual(heat_exact_sol(space,n),ranges) for space in spaces]

def heat_residual(sol,ranges):
    eqn = HeatEquation({'vars':[t,x,y,z],'eqn_kwargs':{'rho':1.,'cp':1.,'k':1.},
                        'sol':[sol],'discontinuities':[]})
    return eqn.balance_integrate(ranges)[0]
    
def test_exact():
    import matplotlib
    import matplotlib.pyplot as plt
    nx = 25
    tests = zip(*heat_exact_tests([[t,0,1],[x,0,1],[y,0,1],[z,0,1]],nx))
    eqns = tests.pop(0)
    angles = tests.pop(0)
    angles = [angle.evalf() for angle in angles[0]]
    tests = list(tests[0])
    for sol in tests:
        plt.plot(angles,sol)
    plt.show()



if __name__=="__main__":
    nx = 10
    vals = numpy.random.rand(10,8)
    vals[:-1,-1] -= .5
    out = []
    ranges = [[t,0,1],[x,0,1],[y,0,1],[z,0,1]]
    for n in range(nx):
        sol = MASA_solution(*(vals[n,:]),cp=1.,k=1.)
        source = lambdify((t,x,y,z),MASA_source(*vals[n,:],cp=1.,k=1.))
        args = (.01,.01,.01,.01)
        dxes = [1. for arg in args]
        out.append(
            abs(recursive_derivative(
                    lambda x0,x1,x2,x3:
                        HeatEquation(sol).balance_integrate(ranges),
                    args,dxes,order=5)-source(*args)))
        print out
    f = open('random_heat_MASA.dat','w')
    f.write('Ax, At, By, Bt, Cz, Ct, Dt, rho')
    for n in range(nx):
        f.write(str(vals[n,:]))
    f.write('')
    f.write(str(out))
    import pdb;pdb.set_trace()
