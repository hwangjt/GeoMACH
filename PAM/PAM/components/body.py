from __future__ import division
from PAM.components import Component, Property
import numpy, pylab, time
import PAM.PAMlib as PAMlib
import mpl_toolkits.mplot3d.axes3d as p3


class Body(Component):
    """ A component used to model blunt bodies. """

    def __init__(self, nx=1, ny=1, nz=1, bottom=2):
        """ Initialization method
        nx: integer
            Number of surfaces in x direction
        ny: integer
            Number of surfaces in y direction
        nz: integer
            Number of surfaces in z direction
        bottom: integer
            Bottom face of the body
            0: open, C0
            1: open, C1
            2: closed
        """ 

        super(Body,self).__init__() 

        self.ms = []
        self.ms.append(numpy.zeros(nx,int))
        self.ms.append(numpy.zeros(ny,int))
        self.ms.append(numpy.zeros(nz,int))

        self.addFace(-2, 3,-0.5)
        self.addFace(-2,-3, 0.5)
        self.addFace( 2, 1,-0.5)
        self.addFace( 3, 1, 0.5)
        self.addFace(-2, 1, 0.5)
        if bottom==2:
            self.addFace(-3, 1,-0.5)

        self.bottom = bottom

    def setDOFs(self):
        setC1 = self.setC1
        setCornerC1 = self.setCornerC1
        for f in range(len(self.Ks)):
            self.setC1('surf', f, val=True)
        if self.bottom==0:
            setC1('surf', 0, i=-1, u=-1, val=False)
            setC1('edge', 0, i=-1, u=-1, val=True)
            setC1('surf', 1, i=-1, u=-1, val=False)
            setC1('edge', 1, i=-1, u=-1, val=True)
            setC1('surf', 2, i= 0, u= 0, val=False)
            setC1('edge', 2, i= 0, u= 0, val=True)
            setC1('surf', 4, i=-1, u=-1, val=False)
            setC1('edge', 4, i=-1, u=-1, val=True)

    def initializeVariables(self):
        nx = self.Qs[2].shape[1]
        ny = self.Qs[2].shape[0]
        nz = self.Qs[3].shape[0]
        zeros = numpy.zeros
        ones = numpy.ones
        self.variables = {
            'noseL':0.1,
            'tailL':0.1,
            'offset':zeros(3),
            'rotation':zeros(3),
            'pos':zeros((nx,3),order='F'),
            'r':zeros((nx,3),order='F'),
            't1U':zeros(nx),
            't2U':zeros(nx),
            't1L':zeros(nx),
            't2L':zeros(nx),
            'shapeL':zeros((ny,nx),order='F'),
            'shapeT':zeros((nz,nx),order='F'),
            'shapeR':zeros((ny,nx),order='F'),
            'shapeB':zeros((nz,nx),order='F')
            }

    def setSections(self, sections=[], t1U=0, t2U=1, t1L=0, t2L=1):
        Ns = self.Ns
        v = self.variables
        for j in range(Ns[2].shape[1]):
            for i in range(Ns[2].shape[0]):
                val = Ns[2][i,j,3]
                if not val == -1:
                    break
            found = False
            for k in range(len(sections)):
                found = found or (val==sections[k])
            if found or sections==[]:
                v['t1U'][j] = t1U
                v['t2U'][j] = t2U
                v['t1L'][j] = t1L
                v['t2L'][j] = t2L

    def propagateQs(self):
        Ns = self.Ns
        Qs = self.Qs
        for f in range(len(Ns)):
            if f==0 or f==4:
                if f==0:
                    L = self.props['noseL']
                    i0 = 0
                    i1 = 1
                    i2 = 2
                else:
                    L = self.props['tailL']
                    i0 = -1
                    i1 = -2
                    i2 = -3
                p = self.props
                dx = abs(self.props['posx'].data[i2] - self.props['posx'].data[i1])  
                Qs[f][:,:,:] = PAMlib.computecone(self.full, Ns[f].shape[0], Ns[f].shape[1], L, p['posy'].data[i0], p['posy'].data[i1], p['posy'].data[i2], p['ry'].data[i1], p['ry'].data[i2], p['rz'].data[i1], p['rz'].data[i2], dx)
                if f==4:
                    Qs[f][:,:,0] *= -1
                    Qs[f][:,:,0] += self.props['noseL']
                    Qs[f][:,:,0] += self.props['posx'].data[-2] - self.props['posx'].data[1]
                    Qs[f][:,:,0] += self.props['tailL']
                    Qs[f][:,:,:] = Qs[f][:,::-1,:]
                Qs[f][:,:,0] += self.offset[0]
                Qs[f][:,:,1] += self.offset[1]
                Qs[f][:,:,2] += self.offset[2]
            else:
                if f==1 and self.full:
                    t1 = 3/4.0
                    t2 = 1/4.0
                elif f==1 and not self.full:
                    t1 = 1/2.0
                    t2 = 1/4.0
                elif f==2:
                    t1 = 1/4.0
                    t2 = -1/4.0
                elif f==3 and self.full:
                    t1 = 7/4.0
                    t2 = 5/4.0
                elif f==3 and not self.full:
                    t1 = 7/4.0
                    t2 = 6/4.0
                elif f==5:
                    t1 = 5/4.0
                    t2 = 3/4.0
                p = self.props
                Qs[f][:,:,:] = PAMlib.computebodysections(Ns[f].shape[0], Ns[f].shape[1], t1, t2, p['noseL'], p['posx'].data, p['posy'].data, p['rz'].data, p['ry'].data, p['t1U'].data, p['t2U'].data, p['t1L'].data, p['t2L'].data)
                Qs[f][:,:,0] += self.offset[0]
                Qs[f][:,:,1] += self.offset[1]
                Qs[f][:,:,2] += self.offset[2]


if __name__ == '__main__':
    #b = Body(nx=2,ny=2,nz=2,bottom=0)
    b = Body(nx=1,ny=2,nz=1,bottom=0)
    import PUBS
    from mayavi import mlab
    b.oml0 = PUBS.PUBS(b.Ps)
    b.setDOFs()
    b.oml0.updateBsplines()
    b.computems()
    b.initializeDOFmappings()
    b.initializeVariables()
    b.setSections()
    #b.propagateQs()
    #b.updateQs()
    b.oml0.computePoints()
    b.oml0.plot(pylab.figure(),False)
    export = PUBS.PUBSexport(b.oml0)
    name='body'
    export.write2Tec(name)
    export.write2TecC(name+'_C')
    pylab.show()
