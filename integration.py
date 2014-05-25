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
'discs' : 
    List of discontinuities for various levels of integration.
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
        self.has_spawned = False
        
    def spawn(self):
        """
Check for additional discontinuities created by integration boundaries.

Returns
-------
iterable
    List of lists of additional discontinuity objects, sorted by level, 
    innermost to outermost.

Notes
-----


"""
        if self.has_spawned:
            raise self.DiscontinuityError(
                "Discontinuity has already spawned.")
        out = []
        # Check intersection of discontinuity with innermost integration 
        # variable:
        var,varmin,varmax = self.ranges[0]
        
        for ind, level in enumerate(self.sym[:-1]):
            # Identify which variables and boundaries apply here.
            
            # Each level is a list of symbolic expressions for the critical 
            # values of the current integration variable (self.vars[ind]).
            new_discs_for_this_level = []
            var = self.vars[ind]
            varmax = self.ranges[ind][2]
            varmin = self.ranges[ind][1]
            for inda,sol in enumerate(level):
                for lim in [varmax, varmin]:
                    new_discs_for_this_level.append(
                        Discontinuity(self.discs[ind][inda].subs({var:lim}),
                                      ranges[ind+1:]))
                    # This returns all solutions. I need a way to eliminate
                    # the imaginary ones.
            out.append(new_discs_for_this_level)
        self.spawn_output = out
        self.has_spawned = True
        return self.spawn_output

    def _local_extrema(self):
        """
Parameters
----------
'disc' : Sympy expression
    Sympy expression f representing discontinuities as f = 0.
'vars' : iterable
    List of sympy variables of integration, in order.

Returns
-------
iterable
    List of lists of sympy expressions for the points where local extrema are 
    to be found, organized by level.

        """
        disc_normal = [sympy.diff(self.disc,var) for var in self.vars]
        eqns = [self.disc]+disc_normal
        self.discs = []
        self.sym = []
        for inda, var in enumerate(self.vars): 
            try:
                sols = sympy.solve(eqns[:inda+1],self.vars[:inda+1],dict=True)
                self.sym.append([sol.pop(var) for sol in sols])
                self.discs.append([self.disc.subs(sol) for sol in sols])
            except(KeyError):
                # No more discontinuities.
                self.sym.extend([[] for indb in range(len(self.vars)-inda-1)])
                break
        return None


    class DiscontinuityError(Exception):
        pass


class Discontinuities(object):
    """List of discontinuities for iterative integration.

Process a list of symbolic discontinuities into a form that is suited for 
numerical integration using SciPy's nquad routine. Specifically, convert a list
of symbolic discontinuities into a list of functions:
[f0(x1,...xn),f1(x2,...xn),...,fn()] where f returns the points of discontinuity
for a given set of arguments.

Parameters
----------
'symbolic_discs' : iterable
    List of symbolic expressions f such that f = 0 is a surface separating two
    distinct states.
'ranges' : iterable
    List of symbolic ranges of the form ((x0,x0_min,x0_max),...,
    (xn,xn_min,xn_max)) where x0,...,xn are Sympy variables, and x_min, x_max
    are numeric values.
'args' : dict, optional
    Dict of any argument substitutions beyond the integration variables.
'opts' : dict, optional
    Unused at present

Attributes
----------

See Also
--------
scipy.integrate.quad : Numerical integration of functions.
scipy.integrate.nquad : Numerical integration of multivariate functions.

Notes
-----
High-accuracy numerical integration of discontinuous, multivariate functions is 
a challenging problem, even if the location of the discontinuity is known 
a priori. Essentially, the problem requires precise fitting of discontinuities,
which is very similar to the problems of shock-fitting and grid generation.
These are hard problems that have no good solutions, and essentially preclude 
the use of truly multidimensional quadrature routines. Iterative integration
based on Fubini's theorem is simpler to apply, provided that all the 
discontinuous points in the various iterated integrals can be identified. That
is what this class does.

Assume a function f0(x0,...,xn), discontinuous at the surface g0(x0,...,xn) = 0,
that must be integrated over some suitable volume V = dx0*...*dxn. Assume also
that integration is carried out in order, such that x0 is the innermost 
integration loop, and xn is the outermost loop. 

The innermost loop is simple: The critical points x0* are found by solving 
g0(x0,...,xn) = 0 for x0 as a function of x1,...xn. The innermost loop can then
be rewritten (assuming only one solution):

f1(x1,...,xn) = 
    int(f0,{dx0,x0L,x0R}) = int(f0,{dx0,x0L,x0*}) + int(f0,{dx0,x0*,x0R})

This type of domain splitting for one-dimensional integrals is provided as part
of the Quadpack library, and therefore also as part of SciPy's quad and nquad
routines using the 'points' keyword. 

The remaining loops are handled similarly, but some care must be taken to 
identify the location of critical points for the new functions being integrated.
Consider the next loop, where f1(x1,...,xn) is integrated over dx1. f1 will have
critical values anywhere df0/dx1 is infinite. This corresponds to one of two
conditions. Either the surface defined by g(x0,...,xn) = 0 is parallel to x0,
or else the surface intersects with the boundary of the integration volume V at 
the bounding surfaces x0L or x0R. 

The first of these cases is simple to detect, and critical points x1* can be 
found by solving the equations {g0 = 0, dg0/dx0 = 0} for {x0, x1} to yield x1* 
as a function of x2,...,xn.

The second case, where the surface g0 = 0 intersects with the bounding surfaces 
x0L or x0R, requires the definition of up to two new discontinuities:
g1L, g1R (x1,...,xn) = g0(x0 = (x0L, x0R), x1,...,xn). These new 
discontinuities have all the same characteristics as the original discontinuity 
g0, except that they are of a lower dimension. They may also be parallel to x1,
just as g0 could be parallel to x0, for instance. They may additionally spawn 
their own discontinuities through intersection with the x1 = x1L, x1R 
boundaries, just as they themselves were spawned by g intersecting the 
x0 = x0L, x0R boundaries.

In general, the critical points xm* of the first kind are given by solving 
{g0 = 0, (dg0/dxi) = 0 | i <= m} for xm, just as was done in the innermost 
case. For instance, x2* is given by solving: 
{g0 = 0, dg0/dx0 = 0, dg0/dx1 = 0, dg0/dx2 = 0} for {x0, x1, x2} and taking the
solution(s) for x2.

Critical values xm* of the second kind arise only if a discontinuity 
intersects with the xmL, xmR boundaries, and these spawn up to two new 
discontinuities.

It is important to note that only discontinuous functions can create 
discontinuities for outer integrals. That is, if once an integration results in
no critical values of either kind, then all subsequent integrations will be 
completely smooth. For best performance, it is therefore good practice to 
choose the order of integration such that discontinuities terminate as quickly 
as possible. 

For instance, suppose a function f0(x0,x1,x2) is discontinuous 
across the surface given by g0(x0,x1,x2) = x0**2 + x1**2 - x2 = 0, and suppose 
we wish to integrate f0 over the range ((x0,-1,1),(x1,-1,1),(x2,-1,2)). The 
discontinuity is parallel to both x0 and x1 at the point (0,0,0), and it also 
intersects with both the x0 and x1 boundaries. 

One could choose to let x0 be the innermost loop and x2 the outermost. If such 
a choice is made, then f1 = int(f0,x0) will have two critical values: one 
formed by the solution of {g0 = 0, dg0/dx0 = 0}, and the other g1 formed by the
intersections of g0 with the x0 boundaries. f2 = int(f1,x1) will have four
xcritical values: One formed by the solution of 
{g0 = 0, dg0/dx0 = 0, dg1/dx1 = 0}, one formed by the solution of 
{g1 = 0, dg1/dx1 = 0}, and two formed by the intersections of g0, g1 with the 
x1 boundaries. If f0 is integrated over an asymmetric domain, then the number 
of discontinuities grows further.

If, on the other hand, one were to choose x2 for the innermost loop, then the 
resulting function f1 = int(f0,x2) would be completely smooth, and no further 
treatment of discontinuities would be necessary. This results in large 
performance increases, as the cost of an integration loop doubles with each
critical point that must be accounted for. 


"""
    def __init__(self,symbolic_discs,ranges,args={},opts={}):
        self.symbolic_discs = symbolic_discs
        self.ranges = ranges
        self.args = args
        self.opts = opts
        self.list_of_discs = [
            Discontinuity(disc,self.ranges,self.args,self.opts) 
            for disc in self.symbolic_discs]
        spawns = [disc.spawn() for disc in self.list_of_discs]
        # Re-group so that like levels are together.
        spawns = [[spawns[inda][indb][indc] 
                   for inda in range(len(spawns))
                   for indc in range(len(spawns[0][0]))
                   ]
                   for indb in range(len(spawns[0]))
                  ]
        levels = []
        levels.append(self._disc_list_constructor(self.list_of_discs))
        no_levels = max([len(spawn) for spawn in spawns])
        for level in spawns:
            levels.append(self._disc_list_constructor(level))
        # Okay, it works up till here. Now I need to re-group again.
        out = []
        for indlvl in range(len(levels)):
            out.append([level.pop(0) for level in levels[:indlvl+1]])
        for level in out:
            # flatten these lists
            level[:] = [item for sublist in level for item in sublist]
        # FINALLY!

    def _disc_list_constructor(self,discs):
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
    # I need to come up with a self-test for this module, especially for the
    # Discontinuity class. It needs to be simple enough that I can reason 
    # through it, yet complex enough to exercise the class. I'm going to use
    # one or more test cases to do this. 
    x=sympy.Symbol('x',real=True)
    y=sympy.Symbol('y',real=True)
    z=sympy.Symbol('z',real=True)
    test_discs = [x**2+y**2+z**2-1,x**2+y**2+z**2-.25]
    ranges = ((x,-.25,1.5),(y,-.25,1.5),(z,-.25,1.5))
    # Begin with a unit sphere that intersects the integration boundary.
    test = Discontinuity(test_discs[0],ranges)
    # Check the discontinuities created by local extrema.
    eqns = []
    # x-values for discontinuity as function of y, z. Order of +- may matter.
    eqns.append([-sympy.sqrt(1-y**2-z**2),sympy.sqrt(1-y**2-z**2)])
    # y-values as function of z. 
    eqns.append([-sympy.sqrt(1-z**2),sympy.sqrt(1-z**2)])
    # z-values.
    eqns.append([-sympy.singleton.S.One,sympy.singleton.S.One])
    if not len(test.sym) == len(eqns):
        raise IntegrationError("""Invalid test solution length.""")
    for level_ind in range(len(test.sym)):
        if not len(test.sym[level_ind]) == len(eqns[level_ind]):
            raise IntegrationError("""Invalid test level length.""")
        for sol_ind in range(len(test.sym[level_ind])):
            x_rand,y_rand,z_rand = [random.random() for item in ranges]
            if not sympy.Equivalent(test.sym[level_ind][sol_ind]==
                                    eqns[level_ind][sol_ind]):
                raise IntegrationError("""
Discontinuity does not return matching 
lower-level non-boundary discontinuities!
""")
    # Check that discontinuities that result from an intersection with a 
    # boundary are also generated correctly.
    test2 = test.spawn()
    # I need some check solutions here, I guess.
    test3 = Discontinuities(test_discs,ranges)
    import pdb;pdb.set_trace()
    phi, theta, psi = .1, .7, .25
    t=sympy.Symbol('t',real=True)
    x=sympy.Symbol('x',real=True)
    y=sympy.Symbol('y',real=True)
    z=sympy.Symbol('z',real=True)
    vars = t,x,y,z
    test = [Discontinuity(x**2+y**2+z**2+t**2-1,ranges),
            Discontinuity(x/t,ranges)]
    test[0].spawn()
    import pdb;pdb.set_trace()
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
