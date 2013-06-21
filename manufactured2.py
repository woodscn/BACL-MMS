import numpy
import sympy
from sympy.utilities.lambdify import lambdify
from sympy.core.cache import clear_cache
import scipy.integrate
from scipy.integrate import nquad

#import nquad

class SympyEquation(object):
    def __init__(self):
        pass
    def dot_product(self,a,b):
        return sum(map(lambda x,y:x*y,a,b))
    def vector_diff(self,expr,var_in):
        vars = list(expr)
        for n in range(len(vars)):
            vars[n] = var_in
        return map(lambda u,x:sympy.diff(u,x),expr,vars)


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
    
class IntegralManufacturedSolution(object):
    pass


class DifferentialIMS(IntegralManufacturedSolution):
    pass


class DifferentialHeatEquation(DifferentialIMS):
    pass


class Error(Exception):
    pass


class IntegrableFunction:

    """
    Inputs:
      SymPy expression for a possibly vector-valued function
      SymPy expressions for various discontinuities
      List of SymPy variables of integration in order e.g. (xi,eta,t)
        ((xi,0,1),(eta,1.3,1.5),(t,-1,1))

    Outputs:
      Function that accepts coordinate values and evaluates the input
      Ranges list for nquad
      Args list for nquad
      Opts list for nquad

    """

    def __init__(
        self,sympy_function,sympy_ranges,sympy_discontinuities=(),args={}):
        
        self.int_variables,self.min_ranges,self.max_ranges=zip(*sympy_ranges)
        self.int_variables = list(self.int_variables)
        self.ranges = zip(self.min_ranges,self.max_ranges)
        self.args=args
        self.sympy_variables = self.int_variables
        self.function = Integrand(
            sympy_function,self.sympy_variables,args=self.args)
        
        int_variables_init = list(self.int_variables)
        discs = list(sympy_discontinuities)
        self.opts = []
        for var in self.int_variables:
            points = []
            int_variables_init.remove(var)
            for disc in discs:
                solved = sympy.solve(disc,var)
                if len(solved) > 0:
                    for solve in solved:
                        points.append(DiscFunction(solve,int_variables_init,
                                                   args=self.args))
                    discs.remove(disc)
            if len(points) > 0:
                self.opts.append(OptionsDict(points))
            else:
                self.opts.append({})
    
    def quad_integrate(self):
        return nquad(self.function,self.ranges,opts=self.opts)

    def mc_integrate(self):
        """
        Evaluate the integral of a function using the Monte Carlo technique.

        mc_integrate(f,ranges,calls,args=())

        Monte Carlo integration is a robust, though inaccurate, method of
        integration. mc_integrate evaluates function f ncalls number of 
        times at random points, and uses the sum of these evaluations to 
        compute the integral of function f over the rectangular volume 
        specified by ranges. The error of the integral scales as 
        1/sqrt(ncalls)

        Inputs:
          f - Callable that returns a number. Must have at least as many 
            number arguments as the number of ranges provided.
          ranges - Sequence of ranges over which to evaluate the integral.
          ncalls - Number of function samples to take. 
          args - Any additional arguments required by f

        Example:
          mc_integrate(f(x,y,z),((0,1),(-.5,.5)),10000,args=3) evaluates
          the integral of f(x,y,z) over the range x=(0,1), y=(-.5,.5), 
          at z=3. The function is sampled 10e4 times, and so the error
          can be expected to be around 0.01.
        """
        f_sum = 0
        for n in xrange(ncalls):
            coords_lst = []
            for rng in ranges:
                coords_lst.append(rng[0] + rng[1]*numpy.random.random())
            f_sum += self.function(coords_lst)
        vol = numpy.prod([float(rng[1])-float(rng[0]) for rng in ranges])
        return vol/calls*f_sum

def int_eqn_sum(eqn_obj,sol,t_range,xi_range,eta_range,zeta_range):
    cons,flux1,flux2,flux3,source = eqn_obj(sol)
    cons_int_lst, flux1_int_lst, flux2_int_lst, flux3_int_lst = [],[],[],[]
    for elem in cons:
        left = IntegrableFunction(elem,(xi_range,eta_range,zeta_range),
                                  sol['discontinuities'],
                                  args={t_range[0]:t_range[1]}
                                  ).quad_integrate()
        right = IntegrableFunction(elem,(xi_range,eta_range,zeta_range),
                                   sol['discontinuities'],
                                   args={t_range[0]:t_range[2]}
                                   ).quad_integrate()
#        left = mc_integrate(cons_sample,(xi_range,eta_range,zeta_range),
#                            calls,args=(t_range[0],elem))
#        right = mc_integrate(cons_sample,(xi_range,eta_range,zeta_range),
#                             calls,args=(t_range[1],elem))
        cons_int_lst.append(right[0]-left[0])
#    print "Done with cons"
    for elem in flux1:
        left = IntegrableFunction(elem,(t_range,eta_range,zeta_range),
                                  sol['discontinuities'],
                                  args={xi_range[0]:xi_range[1]}
                                  ).quad_integrate()
        right = IntegrableFunction(elem,(t_range,eta_range,zeta_range),
                                   sol['discontinuities'],
                                   args={xi_range[0]:xi_range[2]}
                                   ).quad_integrate()
#        left = mc_integrate(flx1_sample,(t_range,eta_range,zeta_range),
#                            calls,args=(xi_range[0],elem))
#        right = mc_integrate(flx1_sample,(t_range,eta_range,zeta_range),
#                             calls,args=(xi_range[1],elem))
        flux1_int_lst.append(right[0]-left[0])
        
#    print "Done with flx1"
    for elem in flux2:
        left = IntegrableFunction(elem,(t_range,xi_range,zeta_range),
                                  sol['discontinuities'],
                                  args={eta_range[0]:eta_range[1]}
                                  ).quad_integrate()
        right = IntegrableFunction(elem,(t_range,xi_range,zeta_range),
                                   sol['discontinuities'],
                                   args={eta_range[0]:eta_range[2]}
                                   ).quad_integrate()
#        left = mc_integrate(flx2_sample,(t_range,xi_range,zeta_range),
#                            calls,(eta_range[0],elem))
#        right = mc_integrate(flx2_sample,(t_range,xi_range,zeta_range),
#                             calls,(eta_range[1],elem))
        flux2_int_lst.append(right[0]-left[0])
#    print "Done with flx2"
    for elem in flux3:
        left = IntegrableFunction(elem,(t_range,xi_range,eta_range),
                                  sol['discontinuities'],
                                  args={zeta_range[0]:zeta_range[1]}
                                  ).quad_integrate()
        right = IntegrableFunction(elem,(t_range,xi_range,eta_range),
                                   sol['discontinuities'],
                                   args={zeta_range[0]:zeta_range[2]}
                                   ).quad_integrate()
#        left = mc_integrate(flx3_sample,(t_range,xi_range,eta_range),
#                            calls,(zeta_range[0],elem))
#        right = mc_integrate(flx3_sample,(t_range,xi_range,eta_range),
#                             calls,(zeta_range[1],elem))
        flux3_int_lst.append(right[0]-left[0])
#    print "Done with flx3"
    source = [0 for elem in cons]
    return (numpy.array(cons_int_lst)+numpy.array(flux1_int_lst)+
            numpy.array(flux2_int_lst)+numpy.array(flux3_int_lst)+
            numpy.array(source))


class OptionsDict(object):
    def __init__(self,points):
        self.points = points
    def __call__(self,*args):
        return {"points":[point(*args) for point in self.points]}


class DiscFunction(object):
    def __init__(self,solved,sympy_variables,args={}):
        self.disc_symbolic = solved.subs(args)
        self.other_sympy_vars = list(sympy_variables)
        self.lambdified = lambdify(self.other_sympy_vars,self.disc_symbolic)
        return None
    def __call__(self,*args):
        if len(args) != len(self.other_sympy_vars):
            print 'args = ',args
            print 'expected args = ',self.other_sympy_vars
            import pdb;pdb.set_trace()
            raise Error('Error in DiscFunction call! Invalid args list!')
        return self.lambdified(*args)


class Integrand(object):
    def __init__(self,sympy_function,sympy_variables,args={}):
        self.sympy_function = sympy_function.subs(args)
        self.sympy_variables = sympy_variables
        self.lambdified = lambdify(self.sympy_variables,self.sympy_function)
        return None
    def __call__(self,*args):
        if len(args) != len(self.sympy_variables):
            print 'args = ',args
            print 'sympy_vars = ',self.sympy_variables
            raise Error('invalid argument list given in call to Integrand!')
        return self.lambdified(*args)

if __name__=="__main__":
    import heat_equation
    t=sympy.Symbol('t')
    xi=sympy.Symbol('xi')
    eta=sympy.Symbol('eta')
    zeta=sympy.Symbol('zeta')
    Ax = 1.5
    nx = 100
    dx = 1./nx
    dy, dz, dt = dx, dx, dx
    ny,nz,nt = 1,1,1
    x0,y0,z0,t0 = 0.,0.,0.,0.
    sources = []
    t_i, t_i_plus_1 = t0, t0+dt
    z_i, z_i_plus_1 = z0, z0+dz
    y_i, y_i_plus_1 = y0, y0+dy
    int_src_array = []
    cos_sol = {}
    cos_sol["solution"] = [sympy.cos(Ax*xi)]
    cos_sol["discontinuities"] = []
    x_array = [x0]
    for i in range(nx+2):
        x_i, x_i_plus_1 = x0 + (i-1)*dx, x0 + (i)*dx
        int_src_array.append(
            int_eqn_sum(heat_equation.HeatEquation(),cos_sol,
                        (t,t_i,t_i_plus_1),(xi,x_i,x_i_plus_1),
                        (eta,y_i,y_i_plus_1),(zeta,z_i,z_i_plus_1)))
        x_array.append(x_i_plus_1)
    x_array = x_array[1:-1]
#    int_src_array = int_src_array[1:-1]
    int_src_array = numpy.array([0]+int_src_array).flatten()
    x_array = numpy.array(x_array)
    src_array = numpy.gradient(numpy.cumsum(int_src_array),dx)[1:-1]/(dy*dz*dt)
    cos_sol_array = Ax**2*scipy.cos(Ax*x_array)
#    print src_array
#    print cos_sol_array
#    print abs(src_array - cos_sol_array)#/abs(cos_sol_array)

