import os
import subprocess
import ctypes
import random
import sympy
import numpy
from sympy.utilities.lambdify import lambdify
from sympy.utilities.codegen import codegen
from sympy.core.cache import clear_cache
from scipy.integrate import nquad
H = sympy.special.delta_functions.Heaviside
from functools import partial

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
        self.integrand = Integrand(
            sympy_function,self.sympy_variables,args=self.args)
        self.function = self.integrand.ctypesified
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

        discs = Discontinuities(sympy_discontinuities,self.sympy_ranges,
                                self.args)
        self.opts = []
        for ind, level in enumerate(self.ranges):
            if discs.nquad_disc_functions:
                self.opts.append(OptionsDict(
                        points=discs.nquad_disc_functions[ind]))
            else:
                self.opts.append({})
        return None

    def quad_integrate(self):
        '''
Integration using scipy.integrate
'''
#        import pdb;pdb.set_trace()
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


class Integrand(object):
    def __init__(self,sympy_function,sympy_variables,args={}):
        self.sympy_function = sympy_function.subs(args)
        self.sympy_variables = sympy_variables
        self.lambdified = lambdify(self.sympy_variables,self.sympy_function)
        clear_cache()
        # I need a unique way to identify the integrand libraries I will be 
        # generating. id(self) works sort of, but it may or may not be unique
        # I think. I'm going to try a hash of the underlying sympy function,
        # but I can't allow negative numbers (C gets confused), so I have to 
        # generate a positive integer instead. Stackoverflow had this idea:
        self.unique_id = ctypes.c_size_t(hash(self.sympy_function)).value
        # stackoverflow.com/questions/18766535/...
        #   positive-integer-from-python-hash-function
        filename_prefix = os.path.join(
            'integrand_libs','integrand'+str(self.unique_id) )
        libname_prefix = os.path.join(
            'integrand_libs',''+str(self.unique_id) )
# Enable the use of ctypes objects in nquad, once multivariate ctypes objects
# are appropriate arguments for the QUADPACK library functions.
        self.generated_code = codegen(
            ('integrand',self.sympy_function),'C',filename_prefix,
            argument_sequence=self.sympy_variables,to_files=True)
        args_str = ",".join(["args["+str(ind)+"]" 
                             for ind in range(len(self.sympy_variables))])
        extra_c_code="".join([r"""
double integrand_wrapper(int n, double args[n])
{
   return integrand(""",args_str,""");

}
"""])
        extra_h_code="""

double integrand_wrapper(int n, double args[n]);

"""
        
        f = open(filename_prefix+".c",'a')
        f.write(extra_c_code)
        f.close()
        f = open(filename_prefix+".h",'a')
        f.write(extra_h_code)
        f.close()
        cmd = ("gcc -dynamiclib -O3 -I. "+filename_prefix+".c -o "
               +filename_prefix+".dylib")
        subprocess.call(cmd,shell=True)
        self.ctypeslib = ctypes.CDLL(filename_prefix+'.dylib')
        self.ctypesified = self.ctypeslib.integrand_wrapper
        self.ctypesified.restype = ctypes.c_double
        # Test the reliability of ctypesified function. This should be 
        # disabled eventually.
        self.ctypesified.argtypes = (ctypes.c_int,
                                     len(self.sympy_variables)*ctypes.c_double)
        
        test = []
        for indx in range(0):
            randargs = [random.random() for item in self.sympy_variables]
            temp = (self.lambdified(*randargs)-
                    self.ctypesified(
                    ctypes.c_int(len(self.sympy_variables)),
                    (len(self.sympy_variables)*ctypes.c_double)(*randargs)))
            if temp**2 >.00000001:
                test.append(temp)
        if test:
            raise IntegrationError("Ctypesified and lambdified do not match!")
        self.ctypesified.argtypes = ctypes.c_int, ctypes.c_double
        return None
    def __call__(self,*args):
        if len(args) != len(self.sympy_variables):
            print 'args = ',args
            print 'sympy_vars = ',self.sympy_variables
            raise Error('invalid argument list given in call to Integrand!')
        import pdb;pdb.set_trace()
        out = self.lambdified(*args)
        out1 = self.ctypesified(len(args),tuple(args))
        print out - out1
        import pdb;pdb.set_trace()
#        print (out-out2)**2
#        exit()
        clear_cache()
        return out


class IntegrationError(Exception):
    pass


class OptionsDict(object):
    def __init__(self,points):
        self.points = points
    def __call__(self,*args):
        if self.points:
            out = {"points":self.points(*args)}
        else:
            out = {}
        return out


class DiscontinuityError(Exception):
    pass


class EmptyDiscontinuity(UserWarning):
    pass


class Discontinuity(object):
    def __init__(self,disc,ranges,args={},opts={}):
        self._disc = disc.subs(args)
        self._ranges = ranges
        self._args = args
        self._opts = opts
        # Eventually set self.method from opts
        self._method = "lambdified"
        self.sym_sols = self._solve_points()
        if not self.sym_sols:
            raise EmptyDiscontinuity()
        self.call_args,self._lambda_list = self._lambdify()
        self.children = self._spawn()

    def __call__(self,*args):
        'Return list of points of discontinuity for given arguments.\n\n'
        if self._method == "lambdified":
            return self._lambdified(*args)
        else:
            raise DiscontinuityError("Undefined call method!")
        
    def __eq__(self,other):
        try:
            out = self._key() == other._key()
        except(AttributeError):
            out = False
        return out

    def __ne__(self,other):
        return not self.__eq__(other)
                
    def _key(self):
        return(type(self).__name__,self._disc,self._ranges,
               self._args,self._opts)

    def __hash__(self):
        return hash(self._key())

    def _solve_points(self):
        try:
            sols = sympy.solve(self._disc,self._ranges[0][0])
        except(KeyError):
            # No solutions.
            sols = []
        return sols
    
    def _lambdify(self):
        lambda_list = []
        vars = [range_[0] for range_ in self._ranges[1:]]
        for sym_sol in self.sym_sols:
            lambda_list.append(lambdify(vars,sym_sol))
        self.__call__.__func__.__doc__ += ('Function signature is f('
                                           +','.join([str(var) for var in vars]
                                                     )+')\n')
        clear_cache()
        return vars,lambda_list
    
    def _lambdified(self,*args):
        return [lambda_(*args) for lambda_ in self._lambda_list]

    def _spawn_local_extrema(self):
        sols = sympy.solve(
            sympy.diff(self._disc,self._ranges[0][0]),self._ranges[0][0])
        new_discs = [self._disc.subs({self._ranges[0][0]:sol}) for sol in sols]
        out = []
        for disc in new_discs:
            try:
                out.append(Discontinuity(disc,self._ranges[1:]))
            except(EmptyDiscontinuity):
                continue
        return out
    
    def _spawn_boundary_intersections(self):
        new_discs = [self._disc.subs({self._ranges[0][0]:lim}) 
                     for lim in self._ranges[0][1:]]
        out = []
        for disc in new_discs:
            try:
                out.append(Discontinuity(disc,self._ranges[1:]))
            except(EmptyDiscontinuity):
                continue
        return out

    def _spawn(self):
        if len(self._ranges) > 1:
            out = (self._spawn_local_extrema() + 
                   self._spawn_boundary_intersections())
        else:
            out = []
        return out


class Discontinuities(object):
    def __init__(self,discs,ranges,args={},opts={}):
        self.ranges = ranges
        self.discs = []
        for disc in discs:
            try:
                self.discs.append(Discontinuity(disc,self.ranges,args,opts))
            except(EmptyDiscontinuity):
                continue
        self.leveled_discs = self._level_discs()
        self.nquad_disc_functions = self.nquad_discs_f()

    def _level_discs(self):
        # Organize the discontinuities according to their level of  
        # integration.
        this_level = list(self.discs)
        out = []
        while not empty(this_level):
            out.append(this_level)
            next_level = []
            for disc in this_level:
                next_level.extend(disc.children)
            this_level = next_level
        # Need to eliminate duplicates
        for level in out:
            new_level = list(level)
            for item in level:
                if new_level.count(item)>1:
                    new_level.remove(item)
            level[:] = new_level
        return out

    class NQuadDiscFunction(object):
        def __init__(self,level):
            self.level = level
        def __call__(self,*args):
            out = []
            for disc in self.level:
                try:
                    out.extend(disc(*args))
                except(ValueError):
                    out.extend([])
            return out

    def nquad_discs_f(self):
        out = [self.NQuadDiscFunction(level) for level in self.leveled_discs]
        return out
        

def empty(seq): # Thanks StackOverflow!
    # See: http://stackoverflow.com/questions/1593564/
    # python-how-to-check-if-a-nested-list-is-essentially-empty
    # Accessed 6 Jun 2014
    try:
        return all(map(empty,seq))
    except TypeError:
        return False


if __name__=="__main__":
    import random
    # Test handling of symbolic discontinuities
    x,y,z = [sympy.Symbol(var,real=True) for var in ['x','y','z']]
    ranges = [[x,-.25,1.25],[y,-.25,1.25],[z,-.25,1.25]]
    sym_disc = x**2+y**2+z**2-1
    disc = Discontinuity(sym_disc,ranges)
    eqns = [-sympy.sqrt(1-y**2-z**2),sympy.sqrt(1-y**2-z**2)]
    for eqn in eqns:
        for sym_sol in disc.sym_sols:
            if sympy.Equivalent(eqn-sym_sol,0):
                break
        else:
            raise DiscontinuityError(
                "Discontinuity test returned incorrect symbolic solutions!")
    yrand,zrand = [.5*random.random()-.25 for ind in [0,1]]
    lambda_sol = disc._lambdified(yrand,zrand)
    subs_sol = [sym_sol.subs({y:yrand,z:zrand}) for sym_sol in disc.sym_sols]
    err = [((lambda_sol[ind]-subs_sol[ind])**2)**.5 
           for ind in range(len(subs_sol))]
    if max(err) > 1*10**-13:
        raise DiscontinuityError(
            "Lambdified solution does not match symbolic solution!")
    test_children = [y**2+z**2-1,y**2+z**2-15./16]
    for test in test_children:
        for child in disc.children:
            if sympy.Equivalent(test,child._disc):
                break
        else:
            raise DiscontinuityError(
                "Direct children do not match!")
    test_discs = [[disc._disc for disc in level] for level in 
                  Discontinuities([sym_disc],ranges).leveled_discs]
    sol_discs = [[x**2+y**2+z**2-1],[y**2+z**2-1,y**2+z**2-15./16],
                 [z**2-1,z**2-15./16,z**2-7./8]]
    for inda in range(len(test_discs)):
        if not set(test_discs[ind])==set(sol_discs[ind]):
            raise DiscontinuityError("Levelled discontinuities do not match!")
    test = Discontinuities([sym_disc],ranges)
    test2 = test.nquad_disc_functions
    args_list = [[.5*random.random()-.25 for ind in range(inds)]
                 for inds in [2,1,0]]
    vars_list = [[y,z],[z],[]]
    for ind, level in enumerate(test.leveled_discs):
        subs_points = [sym_sol.subs(dict(zip(vars_list[ind],args_list[ind])))
               for disc in level for sym_sol in disc.sym_sols]
        subs_points.sort()
        lambda_points = test2[ind](*args_list[ind])
        lambda_points.sort()
        err = [((subs_points[ind] - lambda_points[ind])**2)**.5 < 10**(-13)
               for ind in range(len(subs_points))]
        if [item for item in err if not item]:
            raise DiscontinuityError(
                "Lambdified functions do not match symbolic functions!")

    # Test integration of functions
    test_discs = [x**2+y**2-1]
    # This takes too long at present. Only integrate over x, y.
    test4 = IntegrableFunction(H(test_discs[0]),
                               ranges[:-1],test_discs)
    # Check against evaluation in Mathematica.
    if ((test4.integrate()[0]-0.907360054182972)**2)**.5 > 10**(-7):
        raise IntegrationError(
            "Incorrect integration result for 2-D integration!")


    import pdb;pdb.set_trace()
 
