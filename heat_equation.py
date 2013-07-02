import sympy
from sympy.utilities.lambdify import lambdify

from base_equation import SympyEquation


Zero = sympy.singleton.S.Zero

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

def MASA_solution_full(t_in,x_in,y_in,z_in,Ax,At,By,Bt,Cz,Ct,Dt,rho,cp,k):
    t = sympy.Symbol('t')
    x = sympy.Symbol('x')
    y = sympy.Symbol('y')
    z = sympy.Symbol('z')
    heat_sol = {'vars':[t,x,y,z],'eqn_kwargs':{'rho':rho,'cp':cp,'k':k},
                'sol':[sympy.cos(Ax*x+At*t)*sympy.cos(By*y+Bt*t)*
                       sympy.cos(Cz*z+Ct*t)*sympy.cos(Dt*t)],
                'discontinuities':[]}
    eqn = HeatEquation(heat_sol)
    out = [lambdify(heat_sol['vars'],flux) for flux in eqn.fluxes]
    out.append(lambdify(heat_sol['vars'],eqn.source))
    return out
    
def MASA_solution(t,x,y,z):
    kwargs ={
        'Ax':1., 'At':0., 'By':0., 'Bt':0.,
        'Cz':0., 'Ct':0., 'Dt':0.,
        'rho':1.,'cp':1.,'k':1.
        }
    return MASA_solution_full(t,x,y,z,**kwargs)

if __name__=="__main__":
    t = sympy.Symbol('t')
    x = sympy.Symbol('x')
    y = sympy.Symbol('y')
    z = sympy.Symbol('z')
    # These solutions taken from:
    # http://eqworld.ipmnet.ru/en/solutions/lpde/lpde101.pdf
    A, B, C, mu = 1,1,1,1
    a = 1
    space = sympy.cos(sympy.pi*.25)*x-sympy.sin(sympy.pi*.25)*y
    phi, theta, psi = .1,.7,.25
    space = (sympy.cos(theta)*x-sympy.cos(psi)*sympy.sin(theta)*y+
             sympy.sin(theta)*sympy.sin(psi)*z)
    sols = [A*space+B,
            A*(space**2+2*a*t)+B,
            A*(space**3+6*a*t*space)+B,
            A*(space**4+12*a*t*space**2+12*a**2*t**2)+B,
            A*sympy.exp(a*mu**2*t+mu*space)+B,
            A*sympy.exp(a*mu**2*t-mu*space)+B,
            A*sympy.exp(-a*mu**2*t)*sympy.cos(mu*space+B)+C,
            A*sympy.exp(-mu*space)*sympy.cos(mu*space-2*a*mu**2*t+B)+C]
    error = False
    for sol in sols:
        heat_sol = {'vars':[t,x,y,z],'eqn_kwargs':{'rho':1.,'cp':1.,'k':1.},
                    'sol':[sol],
                    'discontinuities':[]}
        eqn = HeatEquation(heat_sol)
        out = eqn.balance_integrate(((t,0,1),(x,0,1),(y,0,1),(z,0,1)))
        if abs(max(out)) > 1e-14:
            error = True
        print "sol = ", sol
        print "balance integral = ", out
    if error:
        print "Insufficient accuracy in exact solution integral testing!"
    print out
