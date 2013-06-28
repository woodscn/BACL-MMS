from heat_equation import HeatEquation

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
                (self.zeta,self.z_i,self.z_i_plus_1))[0] for j in self.y_ind_array]

def test2():
    p = multiprocessing.Pool()
    dx = 1./Test2MultiClass.nx
    dy, dz, dt = dx, dx, dx
    int_src_array = numpy.array(p.map(test2_multi_func,range(Test2MultiClass.nx+2)))
    src_array = numpy.gradient(numpy.gradient(
            numpy.cumsum(numpy.cumsum(int_src_array,axis=0),axis=1)
            ,dx)[0],dy)[1][1:-1,1:-1]/(dz*dt)
    xgrid = 0*src_array
    for i in range(src_array.shape[0]):
        xgrid[i,:] = (1.-0.)/Test2MultiClass.nx*i
    ygrid = 0*src_array
    for j in range(src_array.shape[1]):
        ygrid[:,j] = (1.-0.)/Test2MultiClass.nx*j
    cos_sol_array = scipy.cos(1.5*xgrid)*scipy.cos(.37*ygrid)
    print "got this far"
    return xgrid,ygrid,(abs(src_array-cos_sol_array)/abs(cos_sol_array))

if __name__=="__main__":
    out = test1()
    plt.plot(out[0],out[1])
    plt.title("1-D, steady-state heat equation, relative error")
    plt.show()
#    out = test2()
    import pdb;pdb.set_trace()

