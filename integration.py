"""
Useful functions for the generation of valid integration calls from Sympy 
expressions.

"""
import sympy
from sympy.utilities.lambdify import lambdify
from sympy.utilities.codegen import codegen
from sympy.utilities.autowrap import autowrap
from sympy.core.cache import clear_cache
from scipy.integrate import nquad
H = sympy.special.delta_functions.Heaviside
import numpy
import subprocess
import ctypes
import itertools
from functools import partial
from multiprocessing import Pool

def list_integral(integrands,**kwargs):
    """
Map numeric integration to a list of symbolic integrands.

Numerically integrate a list of symbolic integrands over a single range and 
with a single set of integration options. 
This is a good spot to implement multiprocessing using Python's 
multiprocessing module.

Parameters
----------
integrands : iterable of sympy expressions
    List of integrand expressions
sympy_ranges : iterable 
    List of ranges and Symbols e.g. ((t,tmin,tmax),(x,xmin,xmax))
sympy_discontinuities : iterable of sympy expressions, optional
    List of symbolic discontinuities e.g. (x/t-0.75). The discontinuity is
    assumed to be at disc == 0, where disc is the sympy expression.
args : dict
    Substitution dict for any additional arguments required for evaluation of
    integrands beyond integration variables given in sympy_ranges.
integrator : not yet implemented
    --
opts : not yet implemented
    --

Returns
-------
iterable of floats
    Results of the integration of integrands over sympy_ranges.
"""
    multi_list_integral = partial(_func_integral,**kwargs)
    out = map(multi_list_integral,integrands)
    # # You have to choose where you will apply multiprocessing, either here or
    # # else in the computation of flux and source integrals. At present, this
    # # is done at the flux/source level.
    # pool = Pool()
    # out = pool.map(multi_list_integral,integrands)
    # pool.close()
    # pool.join()
    return out

def _func_integral(integrand,**kwargs):
    """
Necessary wrapper function to enable use of multiprocessing module.
"""
    out = IntegrableFunction(integrand,**kwargs).integrate()[0]
    return out

class IntegrableFunction(object):
    """
Short Summary: (delete this heading)

Extended Summary: (delete this heading)
Does not allow much in the way of integration options, only providing 'points'.

Parameters
----------
'sympy_function' : Sympy expression
    Sympy expression for a possibly vector-valued function.
'sympy_ranges' : iterable
    List of Sympy variable ranges in order e.g. ((x,0,1),(y,0,1)).
'sympy_discontinuities' : Sympy expression, optional
    Sympy expressions for various discontinuities.
'args' : iterable, optional
    Any additional arguments required by sympy_function.
'integrator' : optional
    Specify an integration method. Unused at present.
'opts' : optional
    Specify function options. Unused at present.

Attributes
----------
'int_variables' : iterable
    List of Sympy Symbols representing the variables of integration.
'ranges' : iterable
    List of ranges, of the form ((xmin,xmax),(ymin,ymax),...).
'function' : lambda or ctypes function 
    Function with signature f(*int_variables,*args).
'integrate' : callable
    Abstracted integration function, currently quadrature.
'symbolic_discontinuities' : iterable
    List of symbolic expressions for discontinuities.
"""

    def __init__(self,sympy_function,sympy_ranges,sympy_discontinuities=(),
                 args={},integrator=None,opts=None):
        self.sympy_ranges = sympy_ranges
        self.int_variables,self.min_ranges,self.max_ranges=zip(*sympy_ranges)
        self.int_variables = list(self.int_variables)
        self.ranges = zip(self.min_ranges,self.max_ranges)
        self.args=args
        self.sympy_variables = self.int_variables
        self.function = Integrand(
            sympy_function,self.sympy_variables,args=self.args).lambdified
        self.integrate = self.quad_integrate

        # Unpack sympy_discontinuities into a list of points for nquad.

        # In order to be used by nquad, discontinuities must be put into a
        # form of functions of the integration variables. One-dimensional 
        # integration effectively smooths a discontinuity, provided the path 
        # of integration crosses the discontinuity. Effectively, this means 
        # that any discontinuity will be smoothed by integration over a 
        # particular variable, provided that the function describing the 
        # discontinuity is dependent on that variable. An example may help. 

        # Assume three discontinuities: [x = 0, x*y = 1, y-1 = 0]. The form of 
        # these discontiuities will depend on the order of integration given 
        # to nquad. If the integration is done as int(int(f(x,y),dx),dy), then 
        # nquad will need the discontinuities in the form: 
        # [[lambda y : 0, lambda y : 1/y],[lambda : 1]]. 
        # Conversely, if the order of integration is reversed to 
        # int(int(f(x,y),dy),dx), then the discontinuities must be
        #  [[lambda x : 1/x, lambda x : 1],[lambda : 0]]. 

        # This segment of code unpacks the list of discontinuities into the 
        # correct form based on the order of integration given by ranges.

        int_variables_init = list(self.int_variables)
        discs = list(sympy_discontinuities)
        self.opts = []
        discontinuity_list = disc_list_constructor(
            [Discontinuity(disc,self.sympy_ranges) for disc in discs])
        self.opts = []
        self.symbolic_discontinuities = []
        vars = list(self.sympy_variables)
        for level in discontinuity_list:
            vars.pop(0)
            self.symbolic_discontinuities.append(
                [solved for solved in level])
            self.opts.append(
                OptionsDict([DiscFunction(solved.as_real_imag()[0],vars)
                             for solved in level]))
        return None
    
    def quad_integrate(self):
        '''
Integration using scipy.integrate
'''
        return nquad(self.function,self.ranges,opts=self.opts)

    def mc_integrate(self):
        """
Integration using Monte Carlo

Monte Carlo integration is a robust, though inaccurate, method of integration. 
mc_integrate evaluates function f ncalls number of times at random points, and 
uses the sum of these evaluations to compute the integral of function f over 
the rectangular volume specified by ranges. The error of the integral scales 
as 1/sqrt(ncalls)
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

class Discontinuity(object):
    """
Representation of a solution discontinuity

Parameters
----------
'disc' : Sympy object
    symbolic representation of the discontinuity as a level-set function
'ranges' : iterable
    List of "range-style" lists, of the form [[x, xmin, xmax],...] where "x" 
    is a sympy variable, and xmin & xmax are floats. Order matters, and 
    eventual integration will assume that the 0-th element is the innermost 
    integral, and so on.
'args' : dict
    Any relevant arguments not contained in 'ranges', e.g. {z:1.101,...} where
    the dict keys must be sympy variables.
'opts' : optional
    Specify additional options. Unused at present.

Attributes
----------
'sym' : 
    List of symbolic expressions for discontinuities to be used in various 
    levels of integration, ordered inner-most to outermost.
'func' :
    List of python lambda function representations of 'sym'.

Methods
-------
'nquad_points'
    Return appropriate DiscFunction objects for high-performance integration in
    nquad.
      
    """
    def __init__(self,disc,ranges,args={},opts={}):
        self.disc = disc
        self.ranges = ranges
        self.args = args
        self.opts = opts
        self.vars = zip(*ranges)[0]
        self._local_extrema()
    
    def _local_extrema(self):
        """
Parameters
----------
'disc' : Sympy expression
    Sympy expression f representing discontinuities as f = 0.
'vars' : iterable
    List of sympy variables of integration, in order.

Returns : iterable
    List of lists of sympy expressions for the points where local extrema are 
    to be found.

        """
        disc_normal = [sympy.diff(self.disc,var) for var in self.vars]
        disc_normal_magnitude_squared = sum([item**2 for item in disc_normal])
        sols = []
        # Compute symbolic expressions for discontinuous points
        for inda in range(len(disc_normal)):
            eqns = [disc_normal[ind] for ind in range(len(disc_normal)-1-inda)]
            eqns.append(self.disc)
            try:
                sols.append([sol[self.vars[-1-inda]] for sol in
                             sympy.solve(eqns,self.vars[:len(self.vars)-inda],
                                         dict=True)])
            except(KeyError):
                sols.append([])
        temp = []
        for sol in sols:
            try:
                temp.append(sol)
            except(IndexError):
                # Don't worry about it if 'sol' is empty.
                pass
        temp.reverse()
        # Compute DiscFunction objects. I'm only interested in real values.
        self.sym = temp
#        import pdb;pdb.set_trace()
#        self.disc_fs = [[DiscFunction(sol.as_real_imag()[0],self.vars[ind+1:])
#                         for sol in sols] for sols in temp]
#        self.sym = [[sol for sol in sols] for sols in temp]
#        # Convert symbolic expressions to python lambda functions
#        self.func = [[lambdify(self.vars[ind+1:],point) for point in points] 
#                     for ind,points in enumerate(self.sym)]
        return None

def disc_list_constructor(discs):
    # Discs is a list of Discontinuity objects, each of which has a sym
    # attribute, which is a list of lists of symbolic point functions.
    # I eventually need a function that returns a list of points for each 
    # level of integration, which means I need to concatenate the innermost
    # lists for each discontinuity, leaving me with a list of lists of
    # symbolic functions, [f(x1,...xn),f(x2,...xn),...f(xn)] that each 
    # contains one or more symbolic expressions of points.
    discontinuities = []
    # Assume each disc is the same length.
    for level in range(len(discs[0].sym)): # For each level of integration.
        # List of discs for a given level
        discontinuities.append(list(
                itertools.chain.from_iterable(
                    [disc.sym[level] for disc in discs])))
    return discontinuities


class DiscFunction(object):
    def __init__(self,solved,sympy_variables,args={}):
        self.disc_symbolic = solved.subs(args)
        self.other_sympy_vars = list(sympy_variables)
        self.lambdified = lambdify(self.other_sympy_vars,self.disc_symbolic)
        clear_cache()
        return None
    def __call__(self,*args):
        if len(args) != len(self.other_sympy_vars):
            print 'args = ',args
            print 'expected args = ',self.other_sympy_vars
            raise Error('Error in DiscFunction call! Invalid args list!')
        out = self.lambdified(*args)
        clear_cache()
        return out


class Integrand(object):
    def __init__(self,sympy_function,sympy_variables,args={}):
        self.sympy_function = sympy_function.subs(args)
        self.sympy_variables = sympy_variables
        self.lambdified = lambdify(self.sympy_variables,self.sympy_function)
        clear_cache()
# # Enable the use of ctypes objects in nquad, once multivariate ctypes objects
# # are appropriate arguments for the QUADPACK library functions.
#        self.generated_code = codegen(
#            ('integrand',self.sympy_function),'C','integrand',
#            argument_sequence=self.sympy_variables,to_files=True)
#        cmd = "gcc -dynamiclib -I. integrand.c -o testlib.dylib"
#        subprocess.call(cmd,shell=True)
#        self.ctypeslib = ctypes.CDLL('./testlib.dylib')
#        self.ctypesified = self.ctypeslib.integrand
#        self.ctypesified.argtypes = tuple(
#            [ctypes.c_double for var in self.sympy_variables])
#        self.ctypesified.restype = ctypes.c_double
        return None
    def __call__(self,*args):
        if len(args) != len(self.sympy_variables):
            print 'args = ',args
            print 'sympy_vars = ',self.sympy_variables
            raise Error('invalid argument list given in call to Integrand!')
        out = self.lambdified(*args)
#        out = self.ctypesified(*args)
#        print (out-out2)**2
#        exit()
        clear_cache()
        return out

class IntegrationError(Exception):
    pass

if __name__=="__main__":
    import random
    H = sympy.special.delta_functions.Heaviside
    phi, theta, psi = .1, .7, .25
    t=sympy.Symbol('t',real=True)
    x=sympy.Symbol('x',real=True)
    y=sympy.Symbol('y',real=True)
    z=sympy.Symbol('z',real=True)
    vars = t,x,y,z
    ranges = ((x,-1.5,1.5),(y,-1.5,1.5),(z,-1.5,1.5),(t,-1.5,1.5))
    test = [Discontinuity(x**2+y**2+z**2+t**2-1,ranges),
            Discontinuity(x/t,ranges)]
    test1 = disc_list_constructor(test)
    opts = []
    vars = [x,y,z,t]
    for level in test1:
        vars.pop(0)
        opts.append(OptionsDict([DiscFunction(solved.as_real_imag()[0],vars)
                     for solved in level]))
    vars = [t,x,y,z]
    test_args = [random.random() for ind in range(len(vars))]
    ranges = ((t,.8,1),(x,-1,1),(y,-1,1),(z,-1,1))
#    S = sympy.functions.Abs(x)/t
    S = x/t
#    S = (sympy.cos(theta)*xi-sympy.cos(psi)*sympy.sin(theta)*eta+
#         sympy.sin(theta)*sympy.sin(psi)*zeta)
#    integrands = [(.5+x**3*t**(-2)-x*y*z)+H(.5-S)+H(.25-S)]
    base_integrand = x**2+y**2-z**2+x*y*z*t+sympy.cos(t)+x*y**3
    integrand = sympy.Piecewise((0+base_integrand,S<0),
                                  (1+base_integrand,True))
    # Sympy has a bug that prevents this integral from being evaluated 
    # correctly. This value comes from Mathematica.
    # The Mathematica code for this is:
    # ranges = {{t, 8/10, 1}, {x, -1, 1}, {y, -1, 1}, {z, -1, 1}};
   # S = x/t;
    # eq = x^2 + y^2 - z^2 + x y z t + Cos[t] + x y^3 + 
    # Piecewise[{{0, S < 0}}, 1];
    # N[Integrate[eq, ranges[[1]], ranges[[2]], ranges[[3]], 
    #             ranges[[4]]], 20]
    exact_integral = 2.3262524846003232935 #sympy.N( 
    #        sympy.integrate(
    #            sympy.integrate(
    #                sympy.integrate(
    #                    sympy.integrate(integrands[0],ranges[0])
    #                    ,ranges[1])
    #                ,ranges[2])
    #            ,ranges[3])
    #        )    
    test = {'sympy_ranges':ranges,'sympy_discontinuities':[S]}#,'args':{z:1}}
    if (abs(IntegrableFunction(integrand,**test).integrate()[0]-exact_integral)
        >1e-14):
        raise IntegrationError("Nquad integration fails for lambda function!")
    ## Ctypes functionality is dependent on improvements to nquad that aren't
    ## yet available.
    # print "ctypes"+str(IntegrableFunction(integrands[0],**test).integrate())
    
    # I'm updating the test cases for use with discontinuities. 
    print "Starting new test case"
    S = t**2+x**2+y**2+z**2 - 0.25
    integrand = 1 - H(S)
    ranges = [[t,-1,1],[x,-1,1],[y,-1,1],[z,-1,1]]
    Discontinuity(S,ranges)
    test = {'sympy_ranges':ranges,'sympy_discontinuities':[S]}
    test = IntegrableFunction(integrand,**test)
    import pdb;pdb.set_trace()
    print test.integrate()
