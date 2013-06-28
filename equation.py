import sympy

class SympyEquation(object):
    def __init__(self):
        pass
    def dot_product(self,a,b):
        return sum(map(lambda x,y:x*y,a,b))
    def vector_diff(self,expr,var_in):
        vars = list(expr)
        for n in range(len(vars)):
            vars[n] = var_in
        return map(lambda u,x:sympy.diff(u,x),expr,vars)


