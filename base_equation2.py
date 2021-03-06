import sympy
import numpy as np
from integration import list_integral
from functools import partial
from multiprocessing import Pool

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
#    def flux_integrate(self,flux,area_ranges,point_range,discs=()):
    def flux_integrate(self,flux,ranges,discs):
        area_ranges_list = [[item for item in ranges if item is not ranges[ind]]
                            for ind in len(ranges)]
        point_ranges_list = [ranges[ind] for ind in range(len(ranges))]
        multi_ranges_list = [
            {'area_ranges':area_ranges_list[ind],
             'point_range_min':
                 {point_ranges_list[ind][0]:point_ranges_list[ind][1]}
             'point_range_max':
                 {point_ranges_list[ind][0]:point_ranges_list[ind][2]}}
            for ind in range(len(ranges))]
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
        partial_junk = partial(junk,obj=self,ranges=ranges,discs=discs)
#        pool = Pool()
        out_list = []
        for ind in range(len(ranges)):
            temp = partial_junk(ind)
            print "Temp is done"
            out_list.append(temp)
#        out_list = map(partial_junk,range(len(ranges)))
#        pool.close()
#        pool.join()
        print "done with flux integrals"
        out = out_list + [self.source_integrate(self.source,ranges,discs)]
        print "done with source integral"
        out = sum(out)
        return out
def junk(ind,area_ranges_list,point_ranges_list,discs):
#    return obj.flux_integrate(
#        obj.fluxes[ind],
#        area_ranges=[item for item in ranges if item is not ranges[ind]],
#        point_range=ranges[ind],discs=discs)

if __name__=="__main__":
    pass
