from itertools import product
import random

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

def heat_exact_sol(**kwargs):
    n, theta, phi, A, B, C, mu, a  = (
        kwargs['n'], kwargs['theta'], kwargs['phi'], kwargs['A'], kwargs['B'],
        kwargs['C'], kwargs['mu'], kwargs['a'] )
    theta, phi = 0.*sympy.pi, 0.*sympy.pi
    space = (sympy.cos(theta)*x
             +sympy.sin(theta)*sympy.cos(phi)*y
             +sympy.sin(theta)*sympy.sin(phi)*z
             )
    # These solutions taken from:
    # http://eqworld.ipmnet.ru/en/solutions/lpde/lpde101.pdf
#    A, B, C, mu = 1,1,1,1
#    a = 1
    sols = [A*space+B,
            A*(space**2+2*a*t)+B,
            A*(space**3+6*a*t*space)+B,
            A*(space**4+12*a*t*space**2+12*a**2*t**2)+B,
            A*sympy.exp(a*mu**2*t+mu*space)+B,
            A*sympy.exp(a*mu**2*t-mu*space)+B,
            A*sympy.exp(-a*mu**2*t)*sympy.cos(mu*space+C)+B,
            A*sympy.exp(-mu*space)*sympy.cos(mu*space-2*a*mu**2*t+C)+B]
    out = {'sol':sympy.Matrix([sols[n]]),'discontinuities':[],
           'vars':[t,x,y,z],'eqn_kwargs':{'rho':1.,'cp':1.,'k':a}}
    return out

def heat_exact_tests(ranges,nx):
    angle_list = [range_*(.5*sympy.pi/(nx-1)) for range_ in range(nx)]
    spaces = [sympy.cos(theta)*x
              +sympy.sin(theta)*sympy.cos(phi)*y
              +sympy.sin(theta)*sympy.sin(phi)*z for angle in angle_list]
    return [[heat_exact_sol(xi,n),angle_list,heat_exact_array(spaces,n,ranges)]
            for n in range(8)]
 

def test_exact():
    ranges = (
        (t,0.,1.)
        ,(x,-.5,.5)
        ,(y,-.5,.5)
        ,(z,-.5,.5)
        )
    n_choices = list(range(8))
    theta_choices = [.05*numpy.pi*range_ for range_ in range(2)]
    phi_choices  = [0]
    A_choices = [.001, .01, .1, 1., 10., 100., 1000.]
    B_choices = A_choices
    C_choices = [.2*range_ for range_ in range(6)]
    mu_choices = [.01, .1, 1., 10., 100.]
    a_choices = A_choices

    n,theta,phi,A,B,C,mu,a = [random.choice(choices) for choices in 
                              [n_choices,theta_choices,phi_choices,A_choices,
                               B_choices,C_choices,mu_choices,a_choices]]
    print n,theta/numpy.pi*180,phi/numpy.pi*180,A,B,C,mu,a
    sol = heat_exact_sol(n=n,theta=theta,phi=phi,A=A,B=B,C=C,mu=mu,a=a)
    S_prime = HeatEquation(sol).balance_integrate(ranges)
    print n,theta/numpy.pi*180,phi/numpy.pi*180,A,B,mu,a,S_prime
    f = open('random_heat_exact.dat','w')
    f.write('problem #, theta(deg), phi(deg), A, B, C, mu, a, residual')
    f.write('\n'+str(
            (n,theta/numpy.pi*180,phi/numpy.pi*180,A,B,C,mu,a,S_prime[0])))
    f.close()

#    import matplotlib
#    import matplotlib.pyplot as plt
#    nx = 25
#    tests = zip(*heat_exact_tests([[t,0,1],[x,0,1],[y,0,1],[z,0,1]],nx))
#    eqns = tests.pop(0)
#    angles = tests.pop(0)
#    angles = [angle.evalf() for angle in angles[0]]
#    tests = list(tests[0])
#    for sol in tests:
#        plt.plot(angles,sol)
#    plt.show()



if __name__=="__main__":
    test_exact()
    print "done"
