""" Conversion utility for triangulation datastructures:
   matplotlib.delaunay.Triangulation
   matplotlib.tri.Triangulation
   CGAL.Constrained_Delaunay_triangulation_2
 See the testing section for how to use these functions.
 requirements: recent matplotlib (1.0.1), numpy
 optional: python-CGAL bindings

 Note regarding triangulations with holes
 matplotlib.delaunay can't deal with holes (definitely in _compute_convex_hull,
 and probably the interpolation code, too).  The work around is to move the
 triangulation to CGAL (where we can enforce the original edges), let CGAL
 fill in edges to make the triangulation convex, then move back to matplotlib.
 The downside is that the interpolation will now fill in any holes.
 There is a bit of a hack here to omit the holes, by labelling the extra triangles
 with a nan inside the interpolator:

 t = tri.Triangulation(x,y,triangles=triangles_with_holes)
 t_complete = convert_triangles.tri_nonconvex_to_convex(t)
 d=convert_triangles.tri_to_delaunay(t_complete)
 lin_interp = d.linear_interpolator(vals)
 lin_interp.planes[len(t.triangles):,0] = nan

 now the interpolation will return nan for holes.


 Rusty Holleman - rustychris@gmail.com
 r2615
"""

import numpy as np
import collections

try:
    from matplotlib import tri
except ImportError:
    # pre 1.0.1 matplotlib doesn't have this module
    print "This module relies on recent matplotlib (1.0.1) with tri module"
    raise

from matplotlib import delaunay

try:
    import CGAL
except:
    CGAL=None

def dummy_delaunay():
    """ For triangulations which are going to be overwritten anyway, feed it
    dummy points which are guaranteed to make a safe triangulation.
    """
    dummy_x = np.array([0.0,0.0,1.0])
    dummy_y = np.array([0.0,1.0,1.0])
    d = delaunay.Triangulation(dummy_x,dummy_y)
    return d

def circumcenter(p1,p2,p3):
    """
    vectorized circumcenter calculation.
    note that this doesn't deal well with nearly degenerate triangles.
    each of p1,p2,p3 is a numpy array of shape [...,2]
    """
    ref = p1

    # translate to reduce roundoff
    p1x = p1[...,0] - ref[...,0] # ==0.0
    p1y = p1[...,1] - ref[...,1] # ==0.0
    p2x = p2[...,0] - ref[...,0]
    p2y = p2[...,1] - ref[...,1]
    p3x = p3[...,0] - ref[...,0]
    p3y = p3[...,1] - ref[...,1]

    vc = np.zeros( p1.shape, np.float64)
    
    # taken from TRANSFORMER_gang.f90
    dd=2.0*((p1x-p2x)*(p1y-p3y) -(p1x-p3x)*(p1y-p2y))
    b1=p1x**2+p1y**2-p2x**2-p2y**2
    b2=p1x**2+p1y**2-p3x**2-p3y**2 
    vc[...,0]=(b1*(p1y-p3y)-b2*(p1y-p2y))/dd + ref[...,0]
    vc[...,1]=(b2*(p1x-p2x)-b1*(p1x-p3x))/dd + ref[...,1]
    
    return vc


def signed_area(x,y):
    """ Signed area, CCW positive, for triangles defined by
    x[...,3], y[...,3].
    i.e. signed_area(x[triangle_nodes],y[triangles]) => [area0,area1,area2,area3...]

    This is not numerically robust.  Use CGAL triangulation if there are nearly degnerate
    triangles.
    """
    i = np.array([0,1,2])
    ip1 = np.array([1,2,0])
    return 0.5*(x[...,i]*y[...,ip1] - x[...,ip1]*y[...,i]).sum(axis=-1)


def tri_to_delaunay(mt,check_ccw=True):
    """ convert matplotlib.tri.Triangulation to matplotlib.delaunay.Triangulation.
    check_ccw: triangles for which signed_area is negative will be reversed.  This
      is necessary to make d._compute_convex_hull() work reliably.  Note that the
      numerics here are not perfect, and if you know that the triangles are already
      ordered CCW, set this to False.

    In general data from mt will be referenced, not copied.
    """
    d = dummy_delaunay()
    d.x = mt.x
    d.y = mt.y
    d.old_shape = d.x.shape

    # copy, in case we have to reorder some vertices
    if check_ccw:
        d.triangle_nodes = mt.triangles.copy()
        areas = signed_area( d.x[d.triangle_nodes],d.y[d.triangle_nodes])
        to_flip = areas < 0
        d.triangle_nodes[to_flip,:] = d.triangle_nodes[to_flip,::-1]
    else:
        d.triangle_nodes = mt.triangles
    
    Ntri = d.triangle_nodes.shape[0]

    P = np.array( [d.x,d.y]).T
    d.circumcenters = circumcenter( P[d.triangle_nodes[:,0]],
                                    P[d.triangle_nodes[:,1]],
                                    P[d.triangle_nodes[:,2]] )

    ab_to_faces = collections.defaultdict(list)

    # build mapping of ordered edge endpoints (a,b) => list of triangles
    for i in range(Ntri):
        # add items so we can look up faces later
        for a,b in [(0,1),(1,2),(2,0)]:
            aa = d.triangle_nodes[i,a]
            bb = d.triangle_nodes[i,b]
            if aa > bb:
                aa,bb = bb,aa
            ab_to_faces[ (aa,bb) ].append(i)

    # Second pass to get triangle neighbors:
    d.triangle_neighbors = np.zeros( (Ntri,3), np.int32) 

    for i in range(Ntri):
        for a,b in [(0,1),(1,2),(2,0)]:
            j = (a+2)%3 # this is edge j, opposite node j, 
            if d.triangle_nodes[i,a] > d.triangle_nodes[i,b]:
                a,b = b,a
            me_and_them = ab_to_faces[ (d.triangle_nodes[i,a],d.triangle_nodes[i,b]) ]

            if len(me_and_them) == 1:
                d.triangle_neighbors[i,j] = -1
            elif me_and_them[0] == i:
                d.triangle_neighbors[i,j] = me_and_them[1]
            else:
                d.triangle_neighbors[i,j] = me_and_them[0]

    d.edge_db = mt.edges
    d.hull = d._compute_convex_hull()
    return d

def tri_to_cgal(mt):
    """ Create a constrained Delaunay triangulation, given the matplotlib.tri.Triangulation
    object mt.

    Notes:
      1. edges are *not* preserved.  If mt is not Delaunay, then the topology will be different.
      2. vertices in the CGAL triangulation have .info() set to the 0-based index of the node.

    Adding constrained edges afterwards:
      cgal_t = tri_to_cgal(...)
      cgal_t.insert_constraint( cgal_t.vh[a], cgal_t.vh[b] )
    """
    Npnt = len(mt.x)
    
    # store vertex handles for future use
    vh = np.zeros( Npnt ,'object')
    DT = CGAL.Constrained_Delaunay_triangulation_2()

    for n in range(Npnt):
        pnt = CGAL.Point_2( mt.x[n], mt.y[n] )
        vh[n] = DT.insert( pnt )
        vh[n].set_info(n)
    DT.vh = vh
    return DT
    
    
def cgal_to_tri(cg,use_info=False):
    """ Convert a CGAL 2-D Delaunay triangulation to a matplotlib.tri.Triangulation
    instance.

    use_info: if False, makes no assumptions that .info() is set on the incoming vertices,
      and will not preserve node numbering.  if True, assumes that vertices in the incoming
      triangulation give the 0-based index of the nodes.

    """
    Npnt = cg.number_of_vertices()

    lookup = {} # (x,y) => index into points
    x = np.zeros( Npnt, np.float64)
    y = np.zeros( Npnt, np.float64)
    triangle_nodes = np.zeros( (cg.number_of_faces(),3), np.int32)
    
    for n,v in enumerate(cg.vertices):
        if use_info:
            n = v.info()
            if n is None:
                raise Exception,"Can't use info() - some vertices are missing data"
        x[n] = v.point().x()
        y[n] = v.point().y()
        if not use_info:
            lookup[ (x[n],y[n]) ] = n

    for i,f in enumerate(cg.faces):
        for k in [0,1,2]:
            pnt = f.vertex(k).point()
            xy = (pnt.x(),pnt.y())
            if use_info:
                n = f.vertex(k).info()
            else:
                n = lookup[xy]
            triangle_nodes[i,k] = n
    return tri.Triangulation(x=x,y=y,triangles=triangle_nodes)
    
def delaunay_to_tri(d):
    """ Return a matplot.tri.Triangulation instance initialized from
    the matplotlib.delaunay.Triangulation object d
    """
    return tri.Triangulation(x=d.x,y=d.y,triangles=d.triangle_nodes)
    
def delaunay_to_cgal(d):
    """ Convert matplotlib.delaunay.Triangulation to CGAL.Constrained_Delaunay_triangulation_2
    """
    return tri_to_cgal( delaunay_to_tri(d))

def cgal_to_delaunay(c):
    """ Convert CGAL.Constrained_Delaunay_triangulation_2 to matplotlib.delaunay.Triangulation """
    return tri_to_delaunay(cgal_to_tri(c))


def tri_to_cgal_nonconvex(mt):
    """
    deal with real-world case, where specified triangulation is non-convex, and has holes.
    the general strategy here is to create a constrained Delaunay triangulation.
    """
    cg = tri_to_cgal(mt)
    for a,b in mt.edges:
        cg.insert_constraint( cg.vh[a], cg.vh[b] )
    return cg

def tri_nonconvex_to_convex(mt):
    # convert to CGAL to fill out the rest of the triangulation
    # as a superset of the original topology
    cg = tri_to_cgal_nonconvex(mt)
    mt_convex = cgal_to_tri(cg,use_info=True)
    Nold = len(mt.triangles)
    Nnew = len(mt_convex.triangles)

    # then reorder to put the original triangles back to the same indices,
    # and move the new, fake triangles to the end of the array

    # the goal is a permutation array such that mt2.triangles[permute,:]
    # will reorder the triangles to match mt.



    # sort each triangle to make it unique
    sorted_new_tris = np.sort( mt_convex.triangles,axis=1 )
    sorted_old_tris = np.sort( mt.triangles,axis=1 )
    mapping = {}
    for i in range(Nnew):
        mapping[ tuple(sorted_new_tris[i]) ] = i
        
    # initialize to -1 so we know whose left over.
    permute = np.zeros( Nnew,np.int32 ) - 1
        
    for i_old in range(Nold):
        i_new = mapping.pop( tuple(sorted_old_tris[i_old] ) )
        # so the new triangle i_new corresponds to the old triangle i_old
        permute[i_old] = i_new

    # everyone who wasn't claimed yet is at the end.
    permute[Nold:] = mapping.values()

    mt_convex.triangles = mt_convex.triangles[permute]
    return mt_convex


####################
# Testing
if __name__ == '__main__':
    Npnts = 50
    pnts = np.random.random((Npnts,2))
    vals = np.random.random(Npnts)

    mt = tri.Triangulation( pnts[:,0],pnts[:,1])

    import pylab
    pylab.figure()

    # tricontour:
    ax1 = pylab.subplot(2,2,1)
    tri.tricontourf(pylab.gca(),mt , vals)
    pylab.axis('equal')
    pylab.title('tri.Triangulation: tricontour')

    # delaunay - linear interpolation
    ax2 =pylab.subplot(2,2,2,sharex=ax1,sharey=ax1)
    d = tri_to_delaunay(mt)
    lin_interp = d.linear_interpolator(vals)
    lin_field = lin_interp[0.0:1.0:100j,0.0:1.0:100j]
    pylab.imshow(lin_field,origin='lower',extent=[0,1,0,1], interpolation='nearest')
    pylab.title('delaunay.Triangulation: linear interp')

    # delaunay - nearest neighbors interpolation
    pylab.subplot(2,2,3,sharex=ax1,sharey=ax1)
    d = tri_to_delaunay(mt)
    nn_interp = d.nn_interpolator(vals)
    nn_field = nn_interp[0.0:1.0:100j,0.0:1.0:100j]
    pylab.imshow(nn_field,origin='lower',extent=[0,1,0,1], interpolation='nearest')
    pylab.title('delaunay.Triangulation: nearest_nbr interp')


    # make the full rounds and show topology
    pylab.subplot(2,2,4,sharex=ax1,sharey=ax1)
    d = tri_to_delaunay(mt)
    if CGAL:
        cg = delaunay_to_cgal(d)
        # just for fun, insert one constrained edge -
        ll = np.argmin(d.x+d.y) ; ur = np.argmax(d.x+d.y)
        cg.insert_constraint( cg.vh[ll], cg.vh[ur] )
        mt3 = cgal_to_tri(cg)
        pylab.title('With and without constrained edge')
    else:
        mt3 = delaunay_to_tri(d)
        pylab.title('Roundtrip to delaunay and back')

    pylab.triplot(mt,lw=2.0,color='b')
    pylab.triplot(mt3,color='r',lw=1.0)
    pylab.axis('equal')
    
    pylab.show()
