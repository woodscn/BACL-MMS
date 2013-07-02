import sympy
from base_equation import SympyEquation

Zero = sympy.singleton.S.Zero

class HeatEquation(SympyEquation):
    def setup(self,**kwargs):
        self.rho = kwargs['rho']
        self.cp = kwargs['cp']
        self.k = kwargs['k']
        self.fluxes = [[self.rho*self.cp*self.sol[0]],
                       [-self.k*sympy.diff(self.sol[0],self.vars_[1])],
                       [-self.k*sympy.diff(self.sol[0],self.vars_[2])],
                       [-self.k*sympy.diff(self.sol[0],self.vars_[3])]]
        self.source = [Zero]
    
if __name__=="__main__":
    Ax = 1. ; At = 0.
    By = 0. ; Bt = 0.
    Cz = 0. ; Ct = 0.
    Dt = 0.
    t = sympy.Symbol('t')
    x = sympy.Symbol('x')
    y = sympy.Symbol('y')
    z = sympy.Symbol('z')
    heat_sol = {'vars':[t,x,y,z],'eqn_kwargs':{'rho':1.,'cp':1.,'k':1.},
                'sol':[sympy.cos(Ax*x+At*t)*sympy.cos(By*y+Bt*t)*
                       sympy.cos(Cz*z+Ct*t)*sympy.cos(Dt*t)],
                'discontinuities':[]}
    eqn = HeatEquation(heat_sol)
    print eqn.balance_integrate(((t,0,1),(x,0,1),(y,0,1),(z,0,1)))
