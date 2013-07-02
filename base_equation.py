import sympy
import numpy as np
from integration import list_integral

class SympyEquation(object):
    def __init__(self,sol):
        self.vars_ = sol['vars']
        self.sol = sol['sol']
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
        return map(lambda u,x:sympy.diff(u,x),expr,vars_)
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
        out = np.zeros(len(self.sol))
        for ind, range_ in enumerate(ranges):
            out = out + self.flux_integrate(
                flux=self.fluxes[ind],
                area_ranges=[item for item in ranges if item is not range_],
                point_range=range_,discs=discs)
        out = out + self.source_integrate(self.source,ranges,discs)
        return out
if __name__=="__main__":
    pass
