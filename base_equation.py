import sympy

class SympyEquation(object):
    def __init__(self,sol):
        self.vars_ = sol['vars']
        self.sol = sol['sol']
        self.setup(**sol['kwargs'])
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

if __name__=="__main__":
    pass
