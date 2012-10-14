from math import sqrt
from cfuncs import tri_area
import cPickle

def integrate_tri(Tri, U):
    """ Integrate values in *U* over triangular grid given by mpl.Triangulation
    instance *Tri*. """
    S = 0.0
    
    for tri in Tri.triangles:
        mean_u = (U[tri[0]] + U[tri[1]] + U[tri[2]]) / 3.0
        area = tri_area((Tri.x[tri[0]], Tri.y[tri[0]]),
                        (Tri.x[tri[1]], Tri.y[tri[1]]),
                        (Tri.x[tri[2]], Tri.y[tri[2]]))
        S += mean_u * area
    return S
    
if __name__ == '__main__':
    with open('testing/tri_data.pkl') as f:
        Tri, U = cPickle.load(f)
    integrate_tri(Tri, U)
    integrate_tri2(Tri, U)
    #integrate_tri3(Tri, U)
    
