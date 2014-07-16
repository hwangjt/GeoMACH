from __future__ import division
import numpy
import pylab
from mpl_toolkits.mplot3d import axes3d
from BSEmodel import BSEmodel


def face(nu, nv, ru, rv, du, dv, d):
    P = numpy.zeros((nu, nv, 3))
    linu = numpy.linspace(-ru, ru, nu)
    linv = numpy.linspace(-rv, rv, nv)
    for i in xrange(nu):
        for j in xrange(nv):
            P[i, j, du] = linu[i]
            P[i, j, dv] = linv[j]
            P[i, j, -du-dv] = d
    return P

def cube(nx, ny, nz, rx, ry, rz):
    return [
        face(nz, ny, rz, ry, 2, 1, -rx),
        face(ny, nz, ry, rz, 1, 2,  rx),
        face(ny, nx, ry, rx, 1, 0, -rz),
        face(nx, ny, rx, ry, 0, 1,  rz),
        face(nx, nz, rx, rz, 0, 2, -ry),
        face(nz, nx, rz, rx, 2, 0,  ry),
        ]

Ps = cube(10, 10, 10, 1, 1, 1)
nsurf = len(Ps)

surf = BSEmodel(Ps)
surf.set_bspline_option('num_pt', 0, 'u', 50)
surf.set_bspline_option('num_pt', 0, 'v', 50)
surf.set_bspline_option('num_pt', 2, 'u', 50)
surf.set_bspline_option('num_pt', 2, 'v', 50)
#for k in xrange(nsurf):
#    surf.set_diff_surf(True, k)
surf.assemble()
surf.print_info()

Cs = cube(4, 4, 4, 1, 1, 1)
for k in xrange(nsurf):
    surf.vec['df_str'](k)[:, :, :] = Cs[k]
surf.apply_jacobian('df', 'd(df)/d(df_str)', 'df_str')
surf.apply_jacobian('cp', 'd(cp)/d(df)', 'df')
surf.apply_jacobian('cp_str', 'd(cp_str)/d(cp)', 'cp')
surf.apply_jacobian('pt_str', 'd(pt_str)/d(cp_str)', 'cp_str')

pts0 = numpy.array([
    [0, 0, 1.2],
    [-1, -1, -1],
    [0.7, -0.8, -0.9],
    ])
surf.compute_projection('test', pts0, ndim=3)
surf.apply_jacobian('test', 'd(test)/d(cp_str)', 'cp_str')
pts = numpy.array(surf.vec['test'].array)
#surf.vec['df'].export_tec_scatter()
surf.vec['pt_str'].export_tec_str()
#surf.vec['pt_str'].export_STL()

if 0:
    Ps = []
    for k in xrange(nsurf):
        Ps.append(surf.vec['pt_str'](k))

    fig = pylab.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(pts0[:, 0], pts0[:, 1], pts0[:, 2])
    ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2])
    for P in Ps:
        ax.plot_wireframe(P[:, :, 0], P[:, :, 1], P[:, :, 2])
    pylab.show()
