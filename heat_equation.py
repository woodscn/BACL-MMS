"""
The linear heat equation, including several solutions
"""
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
    """
Manufactured solution as given in MASA documentation.
"""
    return {'vars':[t,x,y,z],'eqn_kwargs':{'rho':rho,'cp':cp,'k':k},
                'sol':[sympy.cos(Ax*x+At*t)*sympy.cos(By*y+Bt*t)*
                       sympy.cos(Cz*z+Ct*t)*sympy.cos(Dt*t)],
                'discontinuities':[]}

def MASA_source(Ax,At,By,Bt,Cz,Ct,Dt,rho,cp,k):
    """
Analytic manufactured source term as given in MASA documentation.
"""
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
    """
Functional form of 'MASA_source'.

Parameters
----------
(all parameters are keyword arguments)
Ax, At, By, Bt, Cz, Ct, Dt : float
    Parameters used in the computation of 'MASA_source'
rho : float
    mass density of medium.
cp : float
    specific heat at constant pressure of medium.
k : float
    heat conductivity of medium.

Returns
-------
callable
    Returns Python lambda function f(t,x,y,z) to evaluate 'MASA_source' at the
    given point in space and time.
"""
    return lambdify((t,x,y,z),MASA_source(**kwargs))

def heat_exact_sol(**kwargs):
    """
Collection of exact solutions to the heat equation.

These exact solutions are inherently one-dimensional, and are rotated by two
angles in order to obtain solutions that three-dimensional (not aligned with
any three-dimensional coordinate direction).

Parameters
----------
(all parameters are keyword arguments)
n : int in range(8)
    Identifies specific desired exact solution.
theta : float
    Polar angle of solution rotation, measured from x-axis, in radians.
phi : float
    Aximuthal angle of solution rotation, measured from y-axis, in radians.
A, B, C, mu : float
    Exact solution parameters. Can be arbitrary constants.
a : float
    Thermal diffusivity, given by k/(cp*rho)

Returns
-------
dict, containing the fields:
    sol : sympy.Matrix
        Matrix containing symbolic representation of one exact solution.
    discontinuities : list
        Empty list; the linear heat equation does not admit discontinuities.
    vars : list
        List of sympy Symbols: t,x,y,z.
    eqn_kwargs : dict
        kwargs dict that sets 'rho', 'cp', 'k', such that the thermal
        diffusivity 'a' is obtained.
"""
    n, theta, phi, A, B, C, mu, a  = (
        kwargs['n'], kwargs['theta'], kwargs['phi'], kwargs['A'], kwargs['B'],
        kwargs['C'], kwargs['mu'], kwargs['a'] )
#    theta, phi = 0.*sympy.pi, 0.*sympy.pi
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

#def heat_exact_tests(ranges,nx):
#    """
#
#"""
#    angle_list = [range_*(.5*sympy.pi/(nx-1)) for range_ in range(nx)]
#    spaces = [sympy.cos(theta)*x
#              +sympy.sin(theta)*sympy.cos(phi)*y
#              +sympy.sin(theta)*sympy.sin(phi)*z for angle in angle_list]
#    return [[heat_exact_sol(xi,n),angle_list,
#             heat_exact_array(spaces,n,ranges)]
#            for n in range(8)]
 

def test_exact(ntests):
    """
Test the accuracy of numerical integration of the linear heat equation

Compute the weak manufactured source terms from exact solutions for the linear
heat equation, given randomized combinations of the various solution 
parameters. Write the combinations of parameters and the resulting source term
to the file 'random_heat_exact.dat'.

Parameters
----------
ntests : int
    Number of random trials to compute.

Returns
-------
random_heat_exact.dat : ASCII file
    Contains solution parameters and corresponding source term values.
"""
    f = open('random_heat_exact.dat','w')
    f.write('%problem #, theta(deg), phi(deg), A, B, C, mu, a, source')
    ranges = ((t,0.,1.),(x,-.5,.5),(y,-.5,.5),(z,-.5,.5))
    random.seed(100)
    S_prime_list = []
    for indn in range(ntests):
        n_choices = [0,1,2,3,4,5,6,7]
        theta_min, theta_max = 0, numpy.pi*.5
        phi_min, phi_max = 0, numpy.pi*.5
        A_min, A_max = 0.001, 1000
        B_min, B_max = 0.001, 1000
        C_min, C_max = 0, 1
        mu_min, mu_max = 0, 3
        a_min, a_max = 0.001, 10
        n,theta,phi,A,B,C,mu,a = [random.choice(n_choices),
                                  random.random()*(theta_max-theta_min),
                                  random.random()*(phi_max-phi_min),
                                  10**(random.random()*(numpy.log10(A_max)
                                                        -numpy.log10(A_min))
                                       +numpy.log10(A_min)),
                                  10**(random.random()*(numpy.log10(B_max)
                                                        -numpy.log10(B_min))
                                       +numpy.log10(A_min)),
                                  random.random()*(C_max-C_min),
                                  random.random()*(mu_max-mu_min),
                                  10**(random.random()*(numpy.log10(a_max)
                                                        -numpy.log10(a_min))
                                       +numpy.log10(a_min))
                                  ]
        sol = heat_exact_sol(n=n,theta=theta,phi=phi,A=A,B=B,C=C,mu=mu,a=a)
        try:
            S_prime = HeatEquation(sol).balance_integrate(ranges)
            S_prime_list.append(S_prime)
            f.write('\n'+', '.join([str(item) for item in 
                                    (n,theta/numpy.pi*180,phi/numpy.pi*180,
                                     A,B,C,mu,a,S_prime[0])]))
        except(OverflowError):
            print "Overflow Error!"
            print ('n = ',n,'theta = ',theta,'phi = ',phi,
                   'A = ',A,'B = ',B,'C = ',C,'mu = ',mu,'a = ',a)
    f.close()
    return S_prime_list

class HeatEquationError(Exception):
    pass

if __name__=="__main__":
    S_prime_list = test_exact(100)
    S_prime_log = [numpy.log10(numpy.abs(S)) for S in S_prime_list]
    high_error_ratio = (len([S for S in S_prime_log if S > -10])
                        /float(len(S_prime_log)))
    if high_error_ratio > 0.01:
        raise HeatEquationError('Unexpectedly high error in heat equation!')
    # You can load random_heat_exact.dat in Matlab to create a parallel
    # coordinates plot of the error.
    # If you run enough tests, you do find some odd explosions of error for 
    # some cases. These are few (0.33% for 1000 cases with given seed).
    print "All Heat Equation tests passed!"
