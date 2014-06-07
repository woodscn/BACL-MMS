import sympy
from sympy.utilities.lambdify import lambdify
from sympy.core.cache import clear_cache
from functools import partial

class DiscontinuityError(Exception):
    pass


class EmptyDiscontinuity(Exception):
    pass


class Discontinuity(object):
    def __init__(self,disc,ranges,args={},opts={}):
        self._disc = disc
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
        return vars,lambda_list
    
    def _lambdified(self,*args):
        return [lambda_(*args) for lambda_ in self._lambda_list]

    def _spawn_local_extrema(self):
        sols = sympy.solve(sympy.diff(self._disc,self._ranges[0][0]))
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
        self.discs = [Discontinuity(disc,self.ranges,args,opts)
                      for disc in discs]
        self.leveled_discs = self._level_discs()

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

    def nquad_disc_functions(self):
        pass

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
    test.nquad_disc_functions()
    import pdb;pdb.set_trace()
 
