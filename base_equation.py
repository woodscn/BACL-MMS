import sympy
from sympy.utilities.lambdify import lambdify
import numpy as np
from integration import list_integral
from functools import partial
from multiprocessing import Pool

class SympyEquation(object):
    """
Base equation class providing for computation of manufactured source terms.

SympyEquation is an abstract superclass that provides the general routines 
used to compute manufactured source terms given symbolic flux and source terms.

Parameters
----------
sol : dict
    Dictionary containing the entries:
        'vars' - list of sympy Symbols corresponding to integration variables.
        'sol' - sympy Matrix of integrand expressions.
        'discontinuities' - list of sympy expressions corresponding to known 
            solution discontinuities e.g. (x/t-0.75). The discontinuity is
            assumed to be at disc == 0, where disc is the sympy expression.
        'eqn_kwargs' - dict containing any keyword arguments required by the
            specific equation subclass.
    'vars', 'sol', and 'discontinuities' are copied as object attributes.
    'eqn_kwargs' are used as arguments to the setup routine provided by the 
    specific equation subclass.

Attributes
----------
vars : iterable
    List of Sympy Symbols corresponding to integration variables.
sol : sympy.Matrix
    Matrix containing the symbolic expression of the solution to be used. May 
    be a manufactured solution, but may also be an exact solution.
discontinuities : iterable
    List of Sympy expressions corresponding to known solution discontinuities
    e.g. (x/t-0.75). The discontinuity is assumed to be at disc == 0, where
    disc is the sympy expression.
fluxes : list
    List of len('vars') of sympy Matrix objects. Each Matrix is composed of 
    sympy expressions derived from the given solution 'sol' and the flux  
    functions defined by the specific equation subclass and 'eqn_kwargs'.
source : sympy.Matrix
    Sympy Matrix composed of sympy expressions derived from the given solution
    'sol' and the source functions defined by the specific equation subclass
    and 'eqn_kwargs'.
source_diff_symbolic : sympy.Matrix
    symbolic differential source terms, computed by 'balance_diff'


Methods
-------
balance_diff()
    Return symbolic form of strong (differential) manufactured source terms.
balance_lambda_init()
    Return lambdified form of balance_diff.
balance_integrate(ranges,discs=())
    Return weak-form manufactured source terms, numerically integrated over 
    given ranges, with optional discontinuities.
setup(**sol['eqn_kwargs'])
    Subclass-defined initialization routine which computes symbolic 
    expressions of flux and source terms given the solution provided in 'sol'.
    Must set:
        self.fluxes : list of sympy Matrix objects containing sympy 
            expressions.
        self.source : sympy Matrix object containing sympy expressions.
"""
    def __init__(self,sol):
        self.vars_ = sol['vars']
        self.sol = sol['sol']
        self.discontinuities = sol['discontinuities']
        self.setup(**sol['eqn_kwargs'])
    def __call__(self):
        return self.fluxes+[self.source]
    def setup(self):
        pass
    def dot_product(self,a,b):
        return sum(map(lambda x,y:x*y,a,b))
    def vector_diff(self,expr,var_in):
        vars_ = list(expr)
        for n in range(len(vars_)):
            vars_[n] = var_in
#        return sympy.Matrix(map(lambda u,x:sympy.diff(u,x),expr,vars_))
        return sympy.Matrix(map(lambda u,x=var_in:sympy.diff(u,x),expr,vars_))
    def balance_diff(self):
        """
Return symbolic form of strong (differential) manufactured source terms.

Compute strong manufactured source terms based on a differentiable solution,
flux functions, and source functions. The flux functions must be 
differentiable, because symbolic derivatives are taken.

Returns
-------
sympy.Matrix
    Sympy Matrix object containing symbolic manufactured source terms.
"""
        terms = [self.vector_diff(flux,var) 
                 for (flux, var) in zip(self.fluxes,self.vars_)]+[self.source]
        self.source_diff_symbolic = reduce(sympy.Matrix.__add__,terms)
        return self.source_diff_symbolic
    def balance_lambda_init(self):
        """
Lambdify strong symbolic manufactured source terms.

Define python lambda function f(*vars) to evaluate manufactured source terms at
a given point.

Returns
-------
callable
    Python lambda function that evaluates symbolic manufactured source terms at
    a given point (set of 'vars').
"""
        self.balance_diff()
        self.balance_diff_lambda=lambdify(self.vars_,self.source_diff_symbolic)
        return self.balance_diff_lambda
    def flux_integrate(self,flux,area_ranges,point_range,discs=()):
        out = (np.array(list_integral(
                    flux,sympy_ranges=area_ranges,sympy_discontinuities=discs,
                    args={point_range[0]:point_range[2]}))-
               np.array(list_integral(
                    flux,sympy_ranges=area_ranges,sympy_discontinuities=discs,
                    args={point_range[0]:point_range[1]})))
        return out
    def source_integrate(self,source,ranges,discs):
        out = np.array(list_integral(
                source,sympy_ranges=ranges,sympy_discontinuities=discs))
        return out
    def balance_integrate(self,ranges,discs=()):
        """
Evaluate weak (integral) manufactured source terms over a given volume.

Call flux_integrate and source_integrate with the appropriate ranges to 
evaluate the weak (integral) manufactured source terms based on the potentially
discontinuous 'sol' and the flux and source functions. Source terms are 
integrated over the computational volume described by 'ranges'.

Parameters
----------
ranges : iterable
    List of ranges defining the computational volume over which the balance
    integral is to be evaluated e.g. ((t,0,1),(x,-0.5,0.5)) where t,x are sympy
    Symbols. 
discs : iterable
    List of symbolic expressions that describe the location of any known 
    discontinuities in the integrand. The discontinuity is assumed to be
    located where the expressions are equal to zero.

Returns
-------
numpy.array
    Array of floats containing the integrated source terms over 'ranges'.

Notes
-----
This is a promising opportunity for implementing multiprocessing, as the flux-
and source-integrals take some time to compute and are independent of one 
another.
"""
        # Since we're multiprocessing at the list level, we can't
        # also multiprocess here. Pool workers can't spawn their
        # own pools.
        out = np.zeros(len(self.sol))
        wrapper = partial(flux_integrate_wrapper,obj=self,ranges=ranges,
                               discs=self.discontinuities)
#        pool = Pool()
#        out_list = []
#        out_list = pool.map(wrapper,range(len(ranges)))
#        pool.close()
#        pool.join()
        out_list = map(wrapper,range(len(ranges)))
#        print "done with flux integrals"
        out = out_list + [self.source_integrate(self.source,ranges,discs)]
#        print "done with source integral"
        out = sum(out)
        return out
def flux_integrate_wrapper(ind,obj,ranges,discs):
    """
Wrapper for flux_integrate for use in multiprocessing.
"""
    return obj.flux_integrate(
        obj.fluxes[ind],
        area_ranges=[item for item in ranges if item is not ranges[ind]],
        point_range=ranges[ind],discs=discs)

if __name__=="__main__":
    pass
