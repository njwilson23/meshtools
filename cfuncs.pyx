""" Cython utility functions for meshtools """

from math import sqrt

def dot2(tuple a, tuple b):
    cdef float res
    res = a[0]*b[0] + a[1]*b[1]
    return res

def tri_area(tuple a, tuple b, tuple c):
    """ Return the area of a triangle, where (a,b,c) are doubles. """
    cdef tuple ab, ac
    cdef float base, height, p, res

    ab = (b[0]-a[0], b[1]-a[1])
    ac = (c[0]-a[0], c[1]-a[1])
    base = sqrt(ab[0]**2 + ab[1]**2)
    # projection of AC onto AB
    p = dot2(ac, ab) / sqrt(ab[0]**2 + ab[1]**2)
    height = sqrt(ac[0]**2+ac[1]**2 - p**2)
    res = 0.5*base*height
    return res
