import sympy
from sympy.utilities.lambdify import lambdify
import numpy as np
from integration import list_integral
from functools import partial
from multiprocessing import Pool

class SympyEquation(object):
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
        return sympy.Matrix(map(lambda u,x:sympy.diff(u,x),expr,vars_))
    def balance_diff(self):
        terms = [self.vector_diff(flux,var) 
                 for (flux, var) in zip(self.fluxes,self.vars_)]+[self.source]
        self.source_diff_symbolic = reduce(sympy.Matrix.__add__,terms)
        return self.source_diff_symbolic
    def balance_lambda_init(self):
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
        # Since we're multiprocessing at the list level, we can't
        # also multiprocess here. Pool workers can't spawn their
        # own pools.
        out = np.zeros(len(self.sol))
        wrapper = partial(flux_integrate_wrapper,obj=self,ranges=ranges,
                               discs=self.discontinuities)
        pool = Pool()
        out_list = []
        out_list = pool.map(wrapper,range(len(ranges)))
        pool.close()
        pool.join()
        print "done with flux integrals"
        out = out_list + [self.source_integrate(self.source,ranges,discs)]
        print "done with source integral"
        out = sum(out)
        return out
def flux_integrate_wrapper(ind,obj,ranges,discs):
    return obj.flux_integrate(
        obj.fluxes[ind],
        area_ranges=[item for item in ranges if item is not ranges[ind]],
        point_range=ranges[ind],discs=discs)

if __name__=="__main__":
    pass
