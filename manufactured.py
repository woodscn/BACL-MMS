import numpy
import scipy
import sympy
from integration import IntegrableFunction

class Error(Exception):
    pass


#class BaseSolution(object):
#    def __init__(self,*args,**kwargs):
#        self.args = args
#        self.__dict__.update(kwargs)
#
#
class MaterialIntegral(object):
    pass
def flux_volume_integral(eqn_obj,sol,ranges):
eqn
    for flux in fluxes:
        for elem in flux:
            left = IntegrableFunction(elem,
    cons,flux1,flux2,flux3,source = eqn_obj(sol)
    cons_int_lst, flux1_int_lst, flux2_int_lst, flux3_int_lst = [],[],[],[]
    for elem in cons:
        left = IntegrableFunction(elem,(xi_range,eta_range,zeta_range),
                                  sol['discontinuities'],
                                  args={t_range[0]:t_range[1]}
                                  ).integrate()
        right = IntegrableFunction(elem,(xi_range,eta_range,zeta_range),
                                   sol['discontinuities'],
                                   args={t_range[0]:t_range[2]}
                                   ).integrate()
        cons_int_lst.append(right[0]-left[0])
#    print "Done with cons"
    for elem in flux1:
        left = IntegrableFunction(elem,(t_range,eta_range,zeta_range),
                                  sol['discontinuities'],
                                  args={xi_range[0]:xi_range[1]}
                                  ).quad_integrate()
        right = IntegrableFunction(elem,(t_range,eta_range,zeta_range),
                                   sol['discontinuities'],
                                   args={xi_range[0]:xi_range[2]}
                                   ).quad_integrate()
#        left = mc_integrate(flx1_sample,(t_range,eta_range,zeta_range),
#                            calls,args=(xi_range[0],elem))
#        right = mc_integrate(flx1_sample,(t_range,eta_range,zeta_range),
#                             calls,args=(xi_range[1],elem))
        flux1_int_lst.append(right[0]-left[0])
        
#    print "Done with flx1"
    for elem in flux2:
        left = IntegrableFunction(elem,(t_range,xi_range,zeta_range),
                                  sol['discontinuities'],
                                  args={eta_range[0]:eta_range[1]}
                                  ).quad_integrate()
        right = IntegrableFunction(elem,(t_range,xi_range,zeta_range),
                                   sol['discontinuities'],
                                   args={eta_range[0]:eta_range[2]}
                                   ).quad_integrate()
#        left = mc_integrate(flx2_sample,(t_range,xi_range,zeta_range),
#                            calls,(eta_range[0],elem))
#        right = mc_integrate(flx2_sample,(t_range,xi_range,zeta_range),
#                             calls,(eta_range[1],elem))
        flux2_int_lst.append(right[0]-left[0])
#    print "Done with flx2"
    for elem in flux3:
        left = IntegrableFunction(elem,(t_range,xi_range,eta_range),
                                  sol['discontinuities'],
                                  args={zeta_range[0]:zeta_range[1]}
                                  ).quad_integrate()
        right = IntegrableFunction(elem,(t_range,xi_range,eta_range),
                                   sol['discontinuities'],
                                   args={zeta_range[0]:zeta_range[2]}
                                   ).quad_integrate()
#        left = mc_integrate(flx3_sample,(t_range,xi_range,eta_range),
#                            calls,(zeta_range[0],elem))
#        right = mc_integrate(flx3_sample,(t_range,xi_range,eta_range),
#                             calls,(zeta_range[1],elem))
        flux3_int_lst.append(right[0]-left[0])
#    print "Done with flx3"
    source = [0 for elem in cons]
    return (numpy.array(cons_int_lst)+numpy.array(flux1_int_lst)+
            numpy.array(flux2_int_lst)+numpy.array(flux3_int_lst)+
            numpy.array(source))


if __name__=="__main__":
    import heat_equation
    t=sympy.Symbol('t')
    xi=sympy.Symbol('xi')
    eta=sympy.Symbol('eta')
    zeta=sympy.Symbol('zeta')
    Ax = 1.5
    nx = 100
    dx = 1./nx
    dy, dz, dt = dx, dx, dx
    ny,nz,nt = 1,1,1
    x0,y0,z0,t0 = 0.,0.,0.,0.
    sources = []
    t_i, t_i_plus_1 = t0, t0+dt
    z_i, z_i_plus_1 = z0, z0+dz
    y_i, y_i_plus_1 = y0, y0+dy
    int_src_array = []
    cos_sol = {}
    cos_sol["solution"] = [sympy.cos(Ax*xi)]
    cos_sol["discontinuities"] = []
    x_array = [x0]
    for i in range(nx+2):
        x_i, x_i_plus_1 = x0 + (i-1)*dx, x0 + (i)*dx
        int_src_array.append(
            int_eqn_sum(heat_equation.HeatEquation(),cos_sol,
                        (t,t_i,t_i_plus_1),(xi,x_i,x_i_plus_1),
                        (eta,y_i,y_i_plus_1),(zeta,z_i,z_i_plus_1)))
        x_array.append(x_i_plus_1)
    x_array = x_array[1:-1]
#    int_src_array = int_src_array[1:-1]
    int_src_array = numpy.array([0]+int_src_array).flatten()
    x_array = numpy.array(x_array)
    src_array = numpy.gradient(numpy.cumsum(int_src_array),dx)[1:-1]/(dy*dz*dt)
    cos_sol_array = Ax**2*scipy.cos(Ax*x_array)
#    print src_array
#    print cos_sol_array
#    print abs(src_array - cos_sol_array)#/abs(cos_sol_array)

