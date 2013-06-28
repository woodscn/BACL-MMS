import sympy
from sympy.utilities.lambdify import lambdify
from scipy.integrate import nquad
import numpy
from functools import partial
from multiprocessing import Pool

def func_integral(integrand,**kwargs):
    return IntegrableFunction(integrand,**kwargs).integrate()

def list_integral(integrands,**kwargs):
    multi_list_integral = partial(func_integral,**kwargs)
    pool = Pool()
    return pool.map(multi_list_integral,integrands)

def balance_integral(eqn_sol,sympy_ranges,sympy_discs=(),args={},
                     integrator=None):
    assert(len(eqn_sol)-1==len(sympy_ranges))
    fluxes = eqn_sol[0:-1]
    source = eqn_sol[-1]
    flux_with_ind = [(ind,flux) for ind,flux in enumerate(fluxes)]
    pool_flux = partial(flux_integral,ranges=sympy_ranges,discs=sympy_discs,
                        args=args,integrator=integrator)
    fluxes_out = map(pool_flux,flux_with_ind)

class IntegrableFunction(object):

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

    def __init__(self,sympy_function,sympy_ranges,sympy_discontinuities=(),
                 args={},integrator=None):
        self.int_variables,self.min_ranges,self.max_ranges=zip(*sympy_ranges)
        self.int_variables = list(self.int_variables)
        self.ranges = zip(self.min_ranges,self.max_ranges)
        self.args=args
        self.sympy_variables = self.int_variables
        self.function = Integrand(
            sympy_function,self.sympy_variables,args=self.args)
        self.integrate = self.quad_integrate

        # Unpack sympy_discontinuities into a list of points for nquad.

        # In order to be used by nquad, discontinuities must be put into a form 
        # of functions of the integration variables. One-dimensional integration 
        # effectively smooths a discontinuity, provided the path of integration 
        # crosses the discontinuity. Effectively, this means that any 
        # discontinuity will be smoothed by integration over a particular 
        # variable, provided that the function describing the discontinuity is 
        # dependent on that variable. An example may help. 

        # Assume three discontinuities: [x = 0, x*y = 1, y-1 = 0]. The form of 
        # these discontiuities will depend on the order of integration given to 
        # nquad. If the integration is done as int(int(f(x,y),dx),dy), then 
        # nquad will need the discontinuities in the form: [[lambda y : 0, 
        # lambda y : 1/y],[lambda : 1]]. Conversely, if the order of integration 
        # is reversed to int(int(f(x,y),dy),dx), then the discontinuities must 
        # be [[lambda x : 1/x, lambda x : 1],[lambda : 0]]. 

        # This segment of code unpacks the list of discontinuities into the 
        # correct form based on the order of integration given by ranges.

        int_variables_init = list(self.int_variables)
        discs = list(sympy_discontinuities)
        self.opts = []
        for var in self.int_variables:
            points = []
            int_variables_init.remove(var)
            for disc in discs:
                solved = sympy.solve(disc,var)
                if len(solved) > 0: # does disc depend on var?
                    for solve in solved: # add disc to points.
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
    t=sympy.Symbol('t')
    x=sympy.Symbol('x')
    integrands = [t*x-t**2,t**3-x*(t-1)]
    test = {'sympy_ranges':((t,0,1),(x,-1,0))}
    print list_integral(integrands,**test)
    print IntegrableFunction(integrands[0],**test).integrate()
