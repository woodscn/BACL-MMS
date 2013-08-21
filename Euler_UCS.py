import sympy
from sympy.utilities.lambdify import lambdify
H = sympy.special.delta_functions.Heaviside

from base_equation import SympyEquation

Zero = sympy.singleton.S.Zero
t = sympy.Symbol('t')
xi = sympy.Symbol('xi')
eta = sympy.Symbol('eta')
zeta = sympy.Symbol('zeta')
gamma = sympy.Rational(7,5)

#def H(S):
#    if S>0:
#        out = 1.
#    else:
#        if S<0:
#            out = 0.
#        else:
#            out = .5
#    return out

class Euler_UCS(SympyEquation):
    def setup(self,**kwargs):
        (self.pressure,self.density,                # 
         self.vels_x,self.vels_y,self.vels_z,       # u, v, w
         self.dx_dxi,self.dy_dxi,self.dz_dxi,       # A, B, C
         self.dx_deta,self.dy_deta,self.dz_deta,    # L, M, N
         self.dx_dzeta,self.dy_dzeta,self.dz_dzeta, # P, Q, R
         self.dx_dt,self.dy_dt,self.dz_dt,          # U, V, W
         self.x,self.y,self.z) = self.sol
        self.gamma = gamma
        self.jacobian = (
            self.dx_dxi*
            (self.dy_deta*self.dz_dzeta-self.dz_deta*self.dy_dzeta) + 
            self.dy_dxi*
            (self.dz_deta*self.dx_dzeta-self.dx_deta*self.dz_dzeta) +
            self.dz_dxi*
            (self.dx_deta*self.dy_dzeta-self.dy_deta*self.dx_dzeta) )
        self.dxi_dx = (self.dy_deta*self.dz_dzeta-self.dz_deta*self.dy_dzeta
                       )/self.jacobian
        self.dxi_dy = (self.dz_deta*self.dx_dzeta-self.dx_deta*self.dz_dzeta
                       )/self.jacobian
        self.dxi_dz = (self.dx_deta*self.dy_dzeta-self.dy_deta*self.dx_dzeta
                       )/self.jacobian
        self.deta_dx = (self.dz_dxi*self.dy_dzeta-self.dy_dxi*self.dz_dzeta
                        )/self.jacobian
        self.deta_dy = (self.dx_dxi*self.dz_dzeta-self.dz_dxi*self.dx_dzeta
                        )/self.jacobian
        self.deta_dz = (self.dy_dxi*self.dx_dzeta-self.dx_dxi*self.dy_dzeta
                        )/self.jacobian
        self.dzeta_dx = (self.dy_dxi*self.dz_deta-self.dz_dxi*self.dy_deta
                         )/self.jacobian
        self.dzeta_dy = (self.dz_dxi*self.dx_deta-self.dx_dxi*self.dz_deta
                         )/self.jacobian
        self.dzeta_dz = (self.dx_dxi*self.dy_deta-self.dy_dxi*self.dx_deta
                         )/self.jacobian
        self.Dxi_Dt_vec = (
            self.dot_product((self.vels_x-self.dx_dt,self.vels_y-self.dy_dt,
                               self.vels_z-self.dz_dt),
                             (self.dxi_dx,self.dxi_dy,self.dxi_dz)),
            self.dot_product((self.vels_x-self.dx_dt,self.vels_y-self.dy_dt,
                               self.vels_z-self.dz_dt),
                             (self.deta_dx,self.deta_dy,self.deta_dz)),
            self.dot_product((self.vels_x-self.dx_dt,self.vels_y-self.dy_dt,
                               self.vels_z-self.dz_dt),
                             (self.dzeta_dx,self.dzeta_dy,self.dzeta_dz)))
        self.vels_xi_vec = (
            self.dot_product((self.vels_x,self.vels_y,self.vels_z),
                             (self.dxi_dx,self.dxi_dy,self.dxi_dz)),
            self.dot_product((self.vels_x,self.vels_y,self.vels_z),
                             (self.deta_dx,self.deta_dy,self.deta_dz)),
            self.dot_product((self.vels_x,self.vels_y,self.vels_z),
                             (self.dzeta_dx,self.dzeta_dy,self.dzeta_dz)))
        self.energy = (sympy.Rational(1,2)*
                       (self.vels_x**2+self.vels_y**2+self.vels_z**2) + 
                       self.pressure/((self.gamma-1)*self.density))
        self.grad_xi_vec = ((self.dxi_dx,self.dxi_dy,self.dxi_dz),
                            (self.deta_dx,self.deta_dy,self.deta_dz),
                            (self.dzeta_dx,self.dzeta_dy,self.dzeta_dz))
        self.fluxes = [self.cons()]+[self.flux(n) for n in range(3)]
        self.source = self.source_func()
#        out = (self.cons(),self.flux(0),self.flux(1),self.flux(2),self.source_func())
#        return out
    def cons(self):
        return sympy.Matrix([
            self.density*self.jacobian,
            self.density*self.jacobian*self.vels_x,
            self.density*self.jacobian*self.vels_y,
            self.density*self.jacobian*self.vels_z,
            self.density*self.jacobian*self.energy,
            self.dx_dxi,self.dy_dxi,self.dz_dxi,
            self.dx_deta,self.dy_deta,self.dz_deta,
            self.dx_dzeta,self.dy_dzeta,self.dz_dzeta,
            self.dx_dt,self.dy_dt,self.dz_dt,
            self.x,self.y,self.z])
    def flux(self,n):
        out =  [self.jacobian*(self.density*self.Dxi_Dt_vec[n]),
                self.jacobian*(self.density*self.Dxi_Dt_vec[n]*self.vels_x+
                               self.grad_xi_vec[n][0]*self.pressure),
                self.jacobian*(self.density*self.Dxi_Dt_vec[n]*self.vels_y+
                               self.grad_xi_vec[n][1]*self.pressure),
                self.jacobian*(self.density*self.Dxi_Dt_vec[n]*self.vels_z+
                               self.grad_xi_vec[n][2]*self.pressure),
                self.jacobian*(self.density*self.Dxi_Dt_vec[n]*self.energy+
                          self.vels_xi_vec[n]*self.pressure),
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        out[5:8] = [item*sympy.Integer(-1*(1-n)*(1-n/2)) for item in 
                    (self.dx_dt,self.dy_dt,self.dz_dt)]
        out[8:11] = [item*sympy.Integer(-1*(n)*(2-n)) for item in 
                     (self.dx_dt,self.dy_dt,self.dz_dt)]
        out[11:14] = [item*sympy.Integer(-1*(n-1)*(n/2)) for item in 
                      (self.dx_dt,self.dy_dt,self.dz_dt)]
        return sympy.Matrix(out)
    def source_func(self):
        return sympy.Matrix([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                self.dx_dt,self.dy_dt,self.dz_dt])

def unsteady_Euler(case):
    cases = {'simple':simple_case,
             'two_shock':two_shock_case,
#             'one_shock':one_shock_case,
             'normal':normal_case}
    out = {'vars':[t,xi,eta,zeta],'eqn_kwargs':{}}
    out.update(cases[case]())
    return out

def simple_case():
    return {'sol':sympy.Matrix([1,1,1,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0]),
            'discontinuities':[]}

def normal_case():
    theta = 0.
    S = (sympy.cos(theta)*xi+sympy.sin(theta)*eta)/t
    shock_speed = .78931
    p1 = 460.894
    d1 = 5.99924
    umag1 = 19.5975
    M1rel = (umag1-shock_speed)/(gamma*p1/d1)**(.5)
    p2 = ((2.*gamma*M1rel**2-(gamma-1))/(gamma+1))*p1
    d2 = ((gamma+1.)*M1rel**2/((gamma-1.)*M1rel**2+2.))*d1
    umag2 = (1-d1/d2)*shock_speed+umag1*d1/d2
    print "p2,d2,u2 = ",p2,d2,umag2
    u1 = umag1*sympy.cos(theta)
    v1 = umag1*sympy.sin(theta)
    u2 = umag2*sympy.cos(theta)
    v2 = umag2*sympy.sin(theta)
    speeds = [shock_speed]
    states = [sympy.Matrix([p1,d1,u1,v1]),sympy.Matrix([p2,d2,u2,v2])]
    base_state = [0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0]
    out = {'sol':sympy.Matrix(list(states[0]+
                                   H(S-speeds[0])*(states[1]-states[0]))+
                              base_state),
           'discontinuities':[S-speed for speed in speeds]}
    return out

def two_shock_case():
    S = xi/t
    phi, theta, psi = .1, .7, .25
#    S = (sympy.cos(theta)*xi-sympy.cos(psi)*sympy.sin(theta)*eta+
#         sympy.sin(theta)*sympy.sin(psi)*zeta)/t
    speeds = [.789631,8.68975,12.2507]
    states = [sympy.Matrix([460.894,5.99924,19.5975]),
              sympy.Matrix([1691.64,14.2823,8.68975]),
              sympy.Matrix([1691.64,31.0426,8.68975]),
              sympy.Matrix([46.0950,5.99242,-6.19633])]
    base_state = [0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0]
#    out = {'sol':sympy.Matrix(list(H(speeds[0]-S)*states[0]+
#           H(speeds[1]-S)*(1-H(speeds[0]-S))*states[1]+
#           H(speeds[2]-S)*(1-H(speeds[1]-S))*states[2]+
#           (               1-H(speeds[2]-S))*states[3])+base_state),
    out = {'sol':sympy.Matrix(
            list(states[0]+
                 H(S-speeds[0])*(states[1]-states[0])+
                 H(S-speeds[1])*(states[2]-states[1])+
                 H(S-speeds[2])*(states[3]-states[2]))+
            base_state),
           'discontinuities':[S-speed for speed in speeds]}
    return out

def MASA_solution_E():
    kwargs={'x0':1,
            'xx':1,'ax':1,'fx':sympy.sin,
            'xy':1,'ay':1,'fy':sympy.cos,
            'xz':1,'az':1,'fz':sympy.cos,'L':2}
    return {'vars':[t,xi,eta,zeta],'sol':sympy.Matrix(
            [MASA_sol_var(**kwargs) for var in range(5)]+
            [1,0,0,0,1,0,0,0,1,0,0,0,0,0,0]),'discontinuities':[],
            'eqn_kwargs':{}}

def MASA_sol_var(x0,xx,ax,fx,xy,ay,fy,xz,az,fz,L):
    return (x0+
                xx*fx(ax*sympy.pi*xi/L)+
            xy*fy(ay*sympy.pi*eta/L)+
            xz*fz(az*sympy.pi*zeta/L))

if __name__ == "__main__":
    out = MASA_solution_E()
    eqn = Euler_UCS(MASA_solution_E())
    import pdb;pdb.set_trace()
