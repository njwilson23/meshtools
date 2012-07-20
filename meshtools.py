""" Collection of mesh-manipulation and IO functions. """


import sys
import itertools
import datetime
from xml.etree import ElementTree
from xml.dom import minidom
from scipy.spatial import Delaunay
import traceback


#def _parse_coordinates(args):
#    """ Read coordinates in a variety of forms and return a dictionary. """
#
#    # Parse input
#    if not hasattr(args, "__len__"):
#        ndim = len(args[0][0])
#        if ndim >= 2:
#            X = [args[0][i][0] for i in range(len(args[0]))]
#            Y = [args[0][i][1] for i in range(len(args[0]))]
#        else:
#            raise NotImplementedError('ndim={0}'.format(ndim))
#        if ndim >= 3:
#            Z = [args[0][i][2] for i in range(len(args[0]))]
#        else:
#            Z = None
#
#
#    elif len(args) == 1:
#        ndim = len(args[0][0])
#        if ndim >= 2:
#            X = [args[0][i][0] for i in range(len(args[0]))]
#            Y = [args[0][i][1] for i in range(len(args[0]))]
#        else:
#            raise NotImplementedError('ndim={0}'.format(ndim))
#        if ndim >= 3:
#            Z = [args[0][i][2] for i in range(len(args[0]))]
#        else:
#            Z = None
#
#    else
#
#    elif len(args) == 2:
#        ndim = 2
#        (X, Y) = args
#        Z = None
#
#    elif len(args) == 3:
#        ndim = 3
#        (X, Y, Z) = args
#
#    elif len(args) == 4:
#        ndim = 4
#        (X, Y, Z, A) = args
#
#    return {"x":(X,Y,Z), "ndim":ndim}

def read_dolfinxml(fnm):
    """
    Read a Dolfin *.xml mesh file, and return a tuple (V, C) where

    *V*:
        list of vertices in index, coordinate form

    *C*:
        list of cells, with coordinates referencing vertex indices
    """
    XML = ElementTree.ElementTree(file=fnm)
    xmesh = XML.find('mesh')
    dim = int(xmesh.attrib['dim'])
    vertices = []
    cells = []
    II = XML.iter()

    for element in II:
        if element.tag == 'vertex':
            try:
                # Create a iterator of all coord expected based on dim
                c = itertools.chain(
                        [int(element.attrib['index'])],
                        [float(element.attrib['xyz'[a]]) for a in range(dim)])
                # Assemble it and append to vertices
                vertices.append(tuple(c))
            except:
                traceback.print_exc()
                break

        elif element.tag in ('triangle', 'tetrahedron'):
            try:
                # Create a iterator of all coord expected based on dim
                c = itertools.chain(
                        [int(element.attrib['index'])],
                        [int(element.attrib['v'+str(a)]) for a in range(dim+1)])
                # Assemble it and append to cells
                cells.append(tuple(c))
            except:
                traceback.print_exc()
                break

        else:
            pass

    return (vertices, cells)



def write_dolfinxml(fnm, vertices, cells):
    """
    Write a Dolfin *.xml must file.

    *fnm*:
        filename

    *vertices*:
        tuple containing lists of `X`,`Y`... coordinates

    *cells*:
        tuple contianing list of `n0`,`n1`,`n2`... coordinates
    """

    ndim = len(vertices)
    assert ndim in range(4)

    if ndim == 2:
        celltype = "triangle"
    elif ndim == 3:
        celltype = "tetrahedron"

    # Define basic structure
    xml = ElementTree.Element("dolfin",
            attrib={"xmlns:dolfin":"http://fenics.org/dolfin/"})
    xmesh = ElementTree.Element("mesh",
            attrib={"celltype":celltype, "dim":str(ndim)})
    xvertices = ElementTree.Element("vertices",
            attrib={"size":str(len(vertices[0]))})
    xcells = ElementTree.Element("cells", attrib={"size":str(len(cells[0]))})

    # Add individual vertices
    for i in range(len(vertices[0])):
        X = [vertices[j][i] for j in range(len(vertices))]

        attributes = {"index":str(i), "x":str(X[0]), "y":str(X[1])}
        if ndim >= 3:
            attributes["z"] = str(X[2])

        xvertices.append(ElementTree.Element("vertex", attrib=attributes))

    # Add individual cells

    for i in range(len(cells[0])):
        N = [cells[j][i] for j in range(len(cells))]

        attributes = {"index":str(i), "v0":str(N[0]), "v1":str(N[1]),
                        "v2":str(N[2])}
        if ndim >= 3:
            attributes["v3"] = str(N[3])

        xcells.append(ElementTree.Element(celltype, attrib=attributes))

    xmesh.append(xvertices)
    xmesh.append(xcells)
    xml.append(xmesh)

    xmlstring = ElementTree.tostring(xml)

    # Launder through minidom to get a declaration and pretty-printing
    output = minidom.parseString(xmlstring).toprettyxml(indent="  ")

    with open(fnm, "w") as fout:
        fout.write(output)

    return

def read_elefile(fnm):
    """ Read an *.ele file, returning elements. """
    with open(fnm) as fin:
        lines = fin.readlines()

    goodlines = filter(lambda s: s.lstrip()[0] is not "#",
                filter(lambda s: s.strip() is not "", lines))

    n, deg, nattr = goodlines[0].split()
    elementlist = []
    for line in goodlines[1:]:
        elementlist.append(line.split())

    return n, elementlist

def read_nodefile(fnm):
    """ Read a *.node  file, returning vertices. """
    with open(fnm) as fin:
        lines = fin.readlines()

    goodlines = filter(lambda s: s.lstrip()[0] is not "#",
                filter(lambda s: s.strip() is not "", lines))

    n, dim, nattr, nbounds = goodlines[0].split()
    vertexlist = []
    for line in goodlines[1:]:
        vertexlist.append(line.split())

    return n, vertexlist

def write_nodefile(fnm, *args):
    """
    Write a *.node file for reading by triangle.

    Does not support attributes of boundary markers yet.
    """

    # Parse input
    if len(args) == 1:
        ndim = len(args[0][0])
        if ndim >= 2:
            X = [args[0][i][0] for i in range(len(args[0]))]
            Y = [args[0][i][1] for i in range(len(args[0]))]
        else:
            raise NotImplementedError('ndim={0}'.format(ndim))
        if ndim == 3:
            Z = [args[0][i][2] for i in range(len(args[0]))]

    elif len(args) == 2:
        ndim = 2
        (X, Y) = args

    elif len(args) == 3:
        ndim = 3
        (X, Y, Z) = args

    # Check that coordinates have same length
    assert len(X) == len(Y)
    if ndim > 2:
        assert len(X) == len(Z)

    n = len(X)

    # Write *.node file
    with open(fnm, 'w') as fout:

        fout.write("# Written by meshtools.write_nodefile on {0}\n"
                    .format(datetime.datetime.today()))

        fout.write("{n} {dim} {nattr} {nbounds}\n".format(n=n, dim=ndim,
            nattr=0, nbounds=0))

        for i in range(n):
            fout.write("{i} {x} {y}".format(i=i, x=X[i], y=Y[i]))
            if ndim > 2:
                fout.write(" {z}".format(z=Z[i]))
            fout.write("\n")

    return


def read_meshfile(fnm):
    """ Read *.mesh file.

    Based on state:

        0 = read 'Dimension'
        1 = read dimension
        2 = read 'Vertices'
        3 = read number of vertices
        4 = read next vertex
        5 = read 'Triangles' or 'Tetrahedra'
        6 = read number of cells
        7 = read next cell
        8 = done

    based on dolfin.meshconvert.mesh2xml

    """

    print "Reading from Medit format (.mesh)"

    metadata = {}
    vertices = []
    cells = []

    # Open files
    ifile = open(fnm, "r")

    # Scan file for cell type
    cell_type = None
    dim = 0
    while 1:

        # Read next line
        line = ifile.readline()
        if not line: break

        # Remove newline
        if line[-1] == "\n":
            line = line[:-1]

        # Read dimension
        if  line == "Dimension" or line == " Dimension":
            line = ifile.readline()
            metadata['num_dims'] = int(line)
            if metadata['num_dims'] == 2:
                metadata['cell_type'] = "triangle"
                dim = 2
            elif metadata['num_dims'] == 3:
                metadata['cell_type'] = "tetrahedron"
                dim = 3
            break

    # Check that we got the cell type
    if 'cell_type' not in metadata.keys():
        sys.stderr.write("no cell type found\n")
        sys.exit(1)

    # Step to beginning of file
    ifile.seek(0)

    # Current state
    state = 0

    num_vertices_read = 0
    num_cells_read = 0

    while True:

        # Read next line
        line = ifile.readline()
        if not line:
            break

        # Skip comments
        if line[0] == '#':
            continue

        # Remove newline
        if line[-1] == "\n":
            line = line[:-1]

        if state == 0:
            if line == "Dimension" or line == " Dimension":
                state += 1

        elif state == 1:
            #metadata['num_dims'] = int(line)
            state +=1

        elif state == 2:
            if line == "Vertices" or line == " Vertices":
                state += 1

        elif state == 3:
            metadata['num_vertices'] = int(line)
            state +=1

        elif state == 4:
            if metadata['num_dims'] == 2:
                (x, y, tmp) = line.split()
                x = float(x)
                y = float(y)
                z = 0.0
            elif metadata['num_dims'] == 3:
                (x, y, z, tmp) = line.split()
                x = float(x)
                y = float(y)
                z = float(z)
            vertices.append((x, y, z))
            num_vertices_read +=1

            if metadata['num_vertices'] == num_vertices_read:
                state += 1

        elif state == 5:
            if (line == "Triangles"  or line == " Triangles") and metadata['num_dims'] == 2:
                state += 1
            if line == "Tetrahedra" and metadata['num_dims'] == 3:
                state += 1

        elif state == 6:
            metadata['num_cells'] = int(line)
            state +=1

        elif state == 7:
            if metadata['num_dims'] == 2:
                (n0, n1, n2, tmp) = line.split()
                n0 = int(n0) - 1
                n1 = int(n1) - 1
                n2 = int(n2) - 1
                cells.append((n0, n1, n2))
            elif metadata['num_dims'] == 3:
                (n0, n1, n2, n3, tmp) = line.split()
                n0 = int(n0) - 1
                n1 = int(n1) - 1
                n2 = int(n2) - 1
                n3 = int(n3) - 1
                cells.append((n0, n1, n2, n3))
            num_cells_read +=1
            if metadata['num_cells'] == num_cells_read:
                state += 1
        elif state == 8:
            break

    # Check that we got all data
    if state == 8:
        print "Conversion done"
    else:
        sys.stderr.write("missing data\n")

    # Close files
    ifile.close()

    return metadata, vertices, cells

def write_meshfile(fnm, vertices, cells):
    """ Write a *.mesh file. Limited to convex geometries for edge counting. """

    ndim = max((len(a) for a in vertices))
    cell_points = max((len(a) for a in cells))
    if cell_points == 3:
        cell_type = "triangle"
    elif cell_points == 4:
        cell_type = "tetrahedron"
    else:
        raise TypeError('cells must have 3 or 4 vertices')
    ncells = len(cells)
    nvertices = len(vertices)

    try:
        fout = open(fnm, 'w')
        fout.write('MeshVersionFormatted 1\n')
        fout.write('\nDimension\n')
        fout.write(str(ndim) + '\n')

        # Vertices
        fout.write('\nVertices\n')
        fout.write(str(nvertices) + '\n')
        for v in vertices:
            fout.write('{0} {1} {2} {3}\n'.format(v[0], v[1], v[2], 0))

        # Cells
        if cell_type == 'triangle' or ndim == 3:
            fout.write('\nTriangles\n')

            if ndim == 3:
                # Calculate the triangles
                triangles = []
                triangles = list(
                            set(
                             reduce(lambda a,b: a+b,
                             map(lambda a:
                              [tuple(sorted(comb)) for comb in itertools.combinations(a, 3)], cells))))

            else:
                triangles = cells
            fout.write(str(len(triangles)) + '\n')
            for tri in triangles:
                fout.write(' {0} {1} {2} {3}\n'.format(tri[0], tri[1], tri[2], 0))

        if cell_type == 'tetrahedron':
            fout.write('\nTetrahedra\n')
            fout.write(str(len(cells)) + '\n')
            for tet in cells:
                fout.write(' {0} {1} {2} {3} {4}\n'.format(tet[0], tet[1], tet[2], tet[3], 0))

        # Corners
        tessel = Delaunay(vertices)
        print len(tessel.convex_hull)
        corners = []
        for facet in tessel.convex_hull:
            map(corners.append, facet)
        corners = list(set(corners))
        corners.sort()
        fout.write('\nCorners\n')
        fout.write(str(len(corners)) + '\n')
        for corner in corners:
            fout.write(str(corner) + '\n')

        # Edges
        edges = []
        cornerset = set(corners)
        for facet in tessel.convex_hull:
            if len(cornerset.symmetric_difference(facet)) > 0:
                edges.append(sorted(list(cornerset.intersection(facet))))
        fout.write('\nEdges\n')
        fout.write(str(len(edges)) + '\n')
        for edge in edges:
            fout.write(' {0} {1} {2}\n'.format(edge[0], edge[1], 0))

        fout.write('\nEnd')

    except:
        traceback.print_exc()

    finally:
        fout.close()

def write_poly(nodes, fnm=None, output=True):
    """ Take a list of nodes, and write a *.poly file describing the convex hull. """
    tessel = Delaunay(array(nodes))

    # Print #nodes, dimension, attr_bool, bound_bool
    s = "{n} {dim} {attr} {bound}\n".format(n=len(nodes), dim=len(nodes[0]), attr=0, bound=0)

    # List nodes
    for i, node in enumerate(nodes):
        s += " " + str(i)
        s += reduce(lambda a,b:a+b, [" " + str(c) for c in node])
        s += "\n"

    # Print #_facets, bound_attr
    s += "{n} {bound}\n".format(n=len(tessel.convex_hull), bound=0)

    # List facets
    for facet in tessel.convex_hull:
        s += " {npoly} {hole} {bound}\n".format(npoly=1, hole=0, bound=0)
        s += "  {ncorners}".format(ncorners=len(facet))
        s += reduce(lambda a,b:a+b, [" " + str(c) for c in facet])
        s += "\n"

    # Hole list
    s += "0\n"

    # Region list
    s += "0\n"

    if fnm is not None:
        with open(fnm, 'w') as f:
            f.write(s)
    if output:
        return s

    return
