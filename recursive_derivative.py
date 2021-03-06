from functools import partial
from scipy.misc import derivative

def recursive_derivative(func,x0_vec,dx_vec=None,n=1,args=(),order=3):
    if dx_vec is None:
        dx_vec = [1. for elem in x0_vec]
    return _RecursiveDerivative(func,x0_vec,dx_vec,n,order).differentiate(*args)

class _RecursiveDerivative(object):
    def __init__(self,func,x0_vec,dx_vec,n,order):
        self.func=func
        self.x0_vec=x0_vec
        self.dx_vec=dx_vec
        self.n=n
        self.order=order
        return None

    def differentiate(self, *args, **kwargs):
        depth = kwargs.pop('depth', 0)
        if kwargs:
            raise ValueError('unexpected kwargs')
        ind = -depth-1
        x0 = self.x0_vec[ind]
        dx = self.dx_vec[ind]
        if depth + 1 == len(self.x0_vec):
            f = self.func
        else:
            f = partial(self.differentiate, depth=depth+1)
        return derivative(f,x0,dx,self.n,args,self.order)

