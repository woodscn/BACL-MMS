import numpy
import scipy
#import scipy.integrate
import sympy
import matplotlib.pyplot as plt
import multiprocessing

import manufactured

class HeatEquation(manufactured.SympyEquation):
    rho = 1
    cp = 1
    k = 1
    t = sympy.Symbol('t')
    xi = sympy.Symbol('xi')
    eta = sympy.Symbol('eta')
    zeta = sympy.Symbol('zeta')
    def __call__(self,sol):
        self.T = sol["solution"][0]
        out = (self.cons(),self.fluxxi(),self.fluxeta(),self.fluxzeta(),
               self.source())
        return out
    def cons(self):
        return [self.rho*self.cp*self.T]
    def fluxxi(self):
        return [-self.k*sympy.diff(self.T,self.xi)]
    def fluxeta(self):
        return [-self.k*sympy.diff(self.T,self.eta)]
    def fluxzeta(self):
        return [-self.k*sympy.diff(self.T,self.zeta)]
    def source(self):
        return [0]
    
def test1():
    t = sympy.Symbol('t')
    xi = sympy.Symbol('xi')
    eta = sympy.Symbol('eta')
    zeta = sympy.Symbol('zeta')
    Ax = 1.5
    nx = 100
    dx = 1./nx
    dy, dz, dt = dx, dx, dx
    ny, nz, nt = 1, 1, 1
    x0, y0, z0, t0 = 0., 0., 0., 0.
    t_i, t_i_plus_1 = t0, t0+dt
    z_i, z_i_plus_1 = z0, z0+dz
    y_i, y_i_plus_1 = y0, y0+dy
    cos_sol = {}
    cos_sol["solution"] = [sympy.cos(Ax*xi)]
    cos_sol["discontinuities"] = []
    x_array = [x0]
    int_src_array = [0]
    for i in range(nx+2):# It's better to avoid the endpoints with gradient
        x_i, x_i_plus_1 = x0 + (i-1)*dx, x0 + (i)*dx
        int_src_array.append(
            manufactured.int_eqn_sum(HeatEquation(),cos_sol,
                        (t,t_i,t_i_plus_1),(xi,x_i,x_i_plus_1),
                        (eta,y_i,y_i_plus_1),(zeta,z_i,z_i_plus_1)))
        x_array.append(x_i_plus_1)
    int_src_array = numpy.array(int_src_array).flatten()
    x_array = numpy.array(x_array)
    src_array = numpy.gradient(numpy.cumsum(int_src_array),dx)/(dy*dz*dt)
    cos_sol_array = Ax**2*scipy.cos(Ax*x_array)
    print "got this far"
    return x_array[1:-1],(abs(src_array-cos_sol_array)/abs(cos_sol_array))[1:-1]

def test2_multi_func(i):
    return Test2MultiClass()(i)

class Test2MultiClass(object):
    t = sympy.Symbol('t')
    xi = sympy.Symbol('xi')
    eta = sympy.Symbol('eta')
    zeta = sympy.Symbol('zeta')
    nx = 10
    def __init__(self):
        self.Ax = 1.5
        self.By = .37
        self.x0, self.y0, self.z0, self.t0 = 0., 0., 0., 0.
        self.xmax = 1.
        self.dx = (self.xmax-self.x0)/self.nx
        self.dy, self.dz, self.dt = self.dx, self.dx, self.dx
        self.ny, self.nz, self.nt = self.nx, 1, 1
        self.t_i, self.t_i_plus_1 = self.t0, self.t0+self.dt
        self.z_i, self.z_i_plus_1 = self.z0, self.z0+self.dz
        self.cos_sol = {}
        self.cos_sol["solution"] = [
            sympy.cos(self.Ax*self.xi)*sympy.cos(self.By*self.eta)]
        self.cos_sol["discontinuities"] = []
        self.y_ind_array = range(self.ny+2)
    def __call__(self,i):
        print "i = ",i
        return [manufactured.int_eqn_sum(
                HeatEquation(),self.cos_sol,
                (self.t,self.t_i,self.t_i_plus_1),
                (self.xi,self.x0+(i-1)*self.dx,self.x0+(i)*self.dx),
                (self.eta,self.y0+(j-1)*self.dy,self.y0+(j)*self.dy),
                (self.zeta,self.z_i,self.z_i_plus_1)) for j in self.y_ind_array]

def test2():
#    t = sympy.Symbol('t')
#    xi = sympy.Symbol('xi')
#    eta = sympy.Symbol('eta')
#    zeta = sympy.Symbol('zeta')
#    Ax = 1.5
#    By = .37
#    nx = 10
#    dx = 1./nx
#    dy, dz, dt = dx, dx, dx
#    ny, nz, nt = nx, 1, 1
#    x0, y0, z0, t0 = 0., 0., 0., 0.
#    t_i, t_i_plus_1 = t0, t0+dt
#    z_i, z_i_plus_1 = z0, z0+dz
#    cos_sol = {}
#    cos_sol["solution"] = [sympy.cos(Ax*xi)*sympy.cos(By*eta)]
#    cos_sol["discontinuities"] = []
#    int_src_array = numpy.zeros([nx+2,ny+2,nz+2,nt+2])
#    x_array = 0*int_src_array
#    y_array = 0*int_src_array
#    z_array = 0*int_src_array
#    t_array = 0*int_src_array
#    
#    for n in range(nt+2):# It's better to avoid the endpoints with gradient
#        for k in range(nz+2):
#            for j in range(ny+2):
#                print "j = ", j, "of ",ny+2
#                y_i, y_i_plus_1 = y0 + (j-1)*dy, y0+(j)*dy
#                for i in range(nx+2):
#                    x_i, x_i_plus_1 = x0 + (i-1)*dx, x0 + (i)*dx
#                    int_src_array[i,j,k,n] = manufactured.int_eqn_sum(
#                        HeatEquation(),cos_sol,
#                        (t,t_i,t_i_plus_1),(xi,x_i,x_i_plus_1),
#                        (eta,y_i,y_i_plus_1),(zeta,z_i,z_i_plus_1))
#                    x_array[i,j,k,n] = x_i_plus_1
#                    y_array[i,j,k,n] = y_i_plus_1
#                    z_array[i,j,k,n] = z_i_plus_1
#                    t_array[i,j,k,n] = t_i_plus_1
    p = multiprocessing.Pool()
    p.map(test2_multi_func,range(Test2MultiClass.nx+2))
    import pdb;pdb.set_trace()
#    int_src_array = numpy.array(int_src_array).flatten()
    x_array = numpy.array(x_array)
    src_array = numpy.gradient(numpy.cumsum(
            int_src_array[:,1:-1,1,1],0),dx,dy)/(dz*dt)
    cos_sol_array = Ax**2*scipy.cos(Ax*x_array)
    print "got this far"
    return x_array[1:-1],(abs(src_array-cos_sol_array)/abs(cos_sol_array))[1:-1]

if __name__=="__main__":
#    out = test1()
#    plt.plot(out[0],out[1])
#    plt.title("1-D, steady-state heat equation, relative error")
#    plt.show()
    out = test2()

